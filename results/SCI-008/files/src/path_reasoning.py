"""
Explainable path reasoning for drug repurposing predictions.
Finds biological paths between drugs and diseases in the KG.
"""

import os
import json
import pandas as pd
import networkx as nx
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")


def load_graph():
    """Load KG as NetworkX graph."""
    df = pd.read_csv(os.path.join(DATA_DIR, "kg_triples.tsv"), sep="\t")
    with open(os.path.join(DATA_DIR, "entity_types.json")) as f:
        entity_types = json.load(f)

    G = nx.DiGraph()
    for _, row in df.iterrows():
        G.add_node(row["head"], entity_type=entity_types.get(row["head"], "Unknown"))
        G.add_node(row["tail"], entity_type=entity_types.get(row["tail"], "Unknown"))
        G.add_edge(row["head"], row["tail"], relation=row["relation"])
        # Add reverse edges for undirected search
        G.add_edge(row["tail"], row["head"], relation=f"inv_{row['relation']}")

    return G, entity_types


def find_paths(G, source, target, max_length=4):
    """Find all simple paths between source and target up to max_length."""
    paths = []
    try:
        for path in nx.all_simple_paths(G, source, target, cutoff=max_length):
            path_with_relations = []
            for i in range(len(path) - 1):
                edge_data = G.get_edge_data(path[i], path[i+1])
                rel = edge_data.get("relation", "unknown")
                path_with_relations.append({
                    "from": path[i],
                    "relation": rel,
                    "to": path[i+1],
                    "from_type": G.nodes[path[i]].get("entity_type", "Unknown"),
                    "to_type": G.nodes[path[i+1]].get("entity_type", "Unknown"),
                })
            paths.append({
                "nodes": path,
                "length": len(path) - 1,
                "edges": path_with_relations,
            })
    except nx.NetworkXError:
        pass

    return paths


def score_path(path_info, entity_types):
    """Score a path based on biological plausibility."""
    score = 1.0

    # Prefer shorter paths
    score *= (1.0 / path_info["length"])

    # Bonus for paths through key biological entities
    key_relations = {
        "drug_targets_gene": 1.5,
        "gene_associated_disease": 1.5,
        "gene_participates_pathway": 1.3,
        "disease_has_phenotype": 1.2,
        "pathway_involves_phenotype": 1.2,
        "gene_interacts_gene": 1.1,
    }

    for edge in path_info["edges"]:
        rel = edge["relation"]
        base_rel = rel.replace("inv_", "")
        if base_rel in key_relations:
            score *= key_relations[base_rel]

    # Penalty for inverse relations
    inv_count = sum(1 for e in path_info["edges"] if e["relation"].startswith("inv_"))
    score *= (0.8 ** inv_count)

    return score


def interpret_path(path_info):
    """Generate biological interpretation of a path."""
    interpretations = []
    for edge in path_info["edges"]:
        rel = edge["relation"]
        fr = edge["from"]
        to = edge["to"]

        rel_descriptions = {
            "drug_targets_gene": f"{fr} targets {to}",
            "gene_associated_disease": f"{fr} is associated with {to}",
            "gene_participates_pathway": f"{fr} participates in {to}",
            "disease_has_phenotype": f"{to} is a phenotype of {fr}",
            "pathway_involves_phenotype": f"{fr} is involved in {to}",
            "gene_interacts_gene": f"{fr} interacts with {to}",
            "drug_treats_disease": f"{fr} treats {to}",
            "drug_interacts_drug": f"{fr} interacts with {to}",
            "drug_inhibits_gene": f"{fr} inhibits {to}",
            "drug_upregulates_gene": f"{fr} upregulates {to}",
        }

        base_rel = rel.replace("inv_", "")
        if rel.startswith("inv_"):
            desc = f"{to} → {fr} (reverse: {base_rel})"
        else:
            desc = rel_descriptions.get(rel, f"{fr} --[{rel}]--> {to}")

        interpretations.append(desc)

    return " → ".join(interpretations)


def analyze_covid_paths(G, entity_types):
    """Analyze paths from drugs to COVID-19."""
    results = []

    drugs = [e for e, t in entity_types.items() if t == "Drug" and e in G.nodes()]
    target = "COVID-19"

    if target not in G.nodes():
        print("COVID-19 not found in graph")
        return pd.DataFrame()

    for drug in drugs:
        paths = find_paths(G, drug, target, max_length=4)
        for path_info in paths[:5]:  # Top 5 paths per drug
            path_score = score_path(path_info, entity_types)
            interpretation = interpret_path(path_info)
            results.append({
                "drug": drug,
                "disease": target,
                "path_length": path_info["length"],
                "path_score": round(path_score, 4),
                "path_nodes": " → ".join(path_info["nodes"]),
                "interpretation": interpretation,
            })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df = results_df.sort_values("path_score", ascending=False)

    return results_df


def compute_path_statistics(path_results):
    """Compute statistics about discovered paths."""
    if len(path_results) == 0:
        return {}

    stats = {
        "total_paths": len(path_results),
        "unique_drugs": path_results["drug"].nunique(),
        "avg_path_length": path_results["path_length"].mean(),
        "avg_path_score": path_results["path_score"].mean(),
        "path_length_distribution": path_results["path_length"].value_counts().to_dict(),
    }

    # Top drugs by max path score
    top_drugs = (path_results.groupby("drug")["path_score"]
                 .max().sort_values(ascending=False).head(10))
    stats["top_drugs_by_score"] = top_drugs.to_dict()

    return stats


def main():
    print("Loading knowledge graph...")
    G, entity_types = load_graph()
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")

    print("\nAnalyzing paths from drugs to COVID-19...")
    path_results = analyze_covid_paths(G, entity_types)

    if len(path_results) > 0:
        path_results.to_csv(os.path.join(RESULTS_DIR, "covid_paths.csv"), index=False)
        print(f"\nFound {len(path_results)} paths")

        stats = compute_path_statistics(path_results)
        with open(os.path.join(RESULTS_DIR, "path_stats.json"), "w") as f:
            json.dump(stats, f, indent=2, default=str)

        print("\nTop 10 paths by score:")
        top10 = path_results.head(10)
        for _, row in top10.iterrows():
            print(f"  [{row['path_score']:.4f}] {row['drug']} → COVID-19: {row['path_nodes']}")
            print(f"          {row['interpretation']}")
    else:
        print("No paths found.")

    return path_results


if __name__ == "__main__":
    main()
