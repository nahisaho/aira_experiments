"""
Explainable Path Reasoning for Drug Repurposing
Finds biological paths connecting drugs to diseases through genes/pathways
"""

import json
import time
import warnings
from pathlib import Path
from collections import defaultdict

import networkx as nx
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")

BASE = Path(__file__).parent.parent
DATA_DIR = BASE / "data"
RESULTS_DIR = BASE / "results"
LOG_FILE = BASE / "logs" / "process-log.jsonl"


def log_event(event_type, details):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": "EXECUTE",
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": "co-scientist-drug-repurposing",
        **details,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


log_event("handoff_started", {"script": "04_path_reasoning.py"})

# Load KG
df_triples = pd.read_csv(DATA_DIR / "kg_triples.tsv", sep="\t")
entity_df = pd.read_csv(DATA_DIR / "kg_entities.csv")
entity_map = dict(zip(entity_df["id"], entity_df["name"]))
type_map = dict(zip(entity_df["id"], entity_df["type"]))

# Build undirected graph for path search (and directed for reasoning)
G_directed = nx.MultiDiGraph()
G_undirected = nx.MultiGraph()

edge_relations = defaultdict(list)

for _, row in df_triples.iterrows():
    G_directed.add_edge(row["head"], row["tail"], relation=row["relation"])
    G_undirected.add_edge(row["head"], row["tail"], relation=row["relation"])
    edge_relations[(row["head"], row["tail"])].append(row["relation"])
    edge_relations[(row["tail"], row["head"])].append(row["relation"])  # undirected


def get_node_label(node_id):
    name = entity_map.get(node_id, node_id)
    ntype = type_map.get(node_id, "?")
    return f"{name} ({ntype})"


def find_paths(source, target, max_hops=4, max_paths=5):
    """Find all simple paths up to max_hops between source and target."""
    paths = []
    try:
        for path in nx.all_simple_paths(G_undirected, source, target, cutoff=max_hops):
            if len(paths) >= max_paths:
                break
            # Annotate with relations
            annotated = []
            for i in range(len(path)):
                node = path[i]
                annotated.append({
                    "node": node,
                    "label": get_node_label(node),
                    "type": type_map.get(node, "?"),
                })
                if i < len(path) - 1:
                    rels = edge_relations.get((path[i], path[i + 1]), [])
                    if not rels:
                        rels = edge_relations.get((path[i + 1], path[i]), [])
                    annotated.append({"relation": rels[0] if rels else "connected"})
            paths.append(annotated)
    except nx.NetworkXNoPath:
        pass
    return paths


def format_path(annotated_path):
    """Convert annotated path to human-readable string."""
    parts = []
    for item in annotated_path:
        if "node" in item:
            parts.append(item["label"])
        else:
            parts.append(f"--[{item['relation']}]-->")
    return " ".join(parts)


def path_to_dict(annotated_path, path_id):
    """Convert annotated path to serializable dict."""
    nodes = [item["label"] for item in annotated_path if "node" in item]
    relations = [item["relation"] for item in annotated_path if "relation" in item]
    return {
        "path_id": path_id,
        "nodes": nodes,
        "relations": relations,
        "path_str": format_path(annotated_path),
        "length": len(nodes) - 1,
    }


# ─────────────────────────────────────────────
# COVID-19 Drug Repurposing Paths
# ─────────────────────────────────────────────
COVID_ID = "MESH:D000086382"

covid_drugs = {
    "DB14443": "Remdesivir",
    "DB00001X": "Baricitinib",
    "DB00002X": "Tocilizumab",
    "DB01234": "Dexamethasone",
    "DB00010X": "Colchicine",
    "DB00005X": "Hydroxychloroquine",
    "DB01076": "Atorvastatin",
    "DB00795": "Sulfasalazine",
}

all_paths_records = []
path_narratives = []

print("=== Explainable Path Reasoning for COVID-19 ===\n")

for drug_id, drug_name in covid_drugs.items():
    if drug_id not in G_undirected.nodes or COVID_ID not in G_undirected.nodes:
        print(f"  {drug_name}: not in graph")
        continue

    paths = find_paths(drug_id, COVID_ID, max_hops=4, max_paths=3)

    print(f"\n{drug_name} → COVID-19 ({len(paths)} paths found):")
    for i, path in enumerate(paths):
        path_dict = path_to_dict(path, f"{drug_id}_p{i+1}")
        all_paths_records.append({
            "drug_id": drug_id,
            "drug_name": drug_name,
            "disease_id": COVID_ID,
            **path_dict,
        })
        print(f"  Path {i+1}: {path_dict['path_str']}")

        # Generate biological narrative
        nodes = path_dict["nodes"]
        rels = path_dict["relations"]
        if len(nodes) >= 3:
            narrative = (
                f"{drug_name} {rels[0]} {nodes[1]}, "
                f"which is {rels[1] if len(rels) > 1 else 'connected to'} {nodes[2] if len(nodes) > 2 else 'COVID-19'}"
            )
            if len(nodes) > 3:
                narrative += f", ultimately affecting COVID-19"
            path_narratives.append({
                "drug": drug_name,
                "narrative": narrative,
                "path_length": path_dict["length"],
            })

# ─────────────────────────────────────────────
# Meta-path analysis
# ─────────────────────────────────────────────
print("\n\n=== Meta-path Statistics ===")

meta_paths = defaultdict(int)
for rec in all_paths_records:
    rels = tuple(rec.get("relations", []))
    if rels:
        meta_paths[rels] += 1

# Top meta-paths
meta_df = pd.DataFrame([
    {"meta_path": " → ".join(k), "count": v}
    for k, v in sorted(meta_paths.items(), key=lambda x: -x[1])
])
if not meta_df.empty:
    print(meta_df.head(10).to_string())

# ─────────────────────────────────────────────
# Path importance scoring
# ─────────────────────────────────────────────
# Score paths by node centrality (hub genes are more important)
centrality = nx.degree_centrality(G_undirected)

for rec in all_paths_records:
    nodes_raw = [
        item["node"] for item in 
        find_paths(rec["drug_id"], rec["disease_id"], max_hops=4, max_paths=1)[0]
        if "node" in item
    ] if find_paths(rec["drug_id"], rec["disease_id"], max_hops=4, max_paths=1) else []
    
    intermediate_centrality = np.mean([
        centrality.get(n, 0)
        for n in nodes_raw[1:-1]
    ]) if nodes_raw[1:-1] else 0
    
    rec["intermediate_centrality"] = round(intermediate_centrality, 4)

# Save results
df_paths = pd.DataFrame(all_paths_records)
df_paths.to_csv(RESULTS_DIR / "drug_disease_paths.csv", index=False)

df_narratives = pd.DataFrame(path_narratives)
if not df_narratives.empty:
    df_narratives.to_csv(RESULTS_DIR / "path_narratives.csv", index=False)
    print("\n\n=== Biological Narratives ===")
    for _, row in df_narratives.iterrows():
        print(f"  [{row['drug']}] {row['narrative']}")

if not meta_df.empty:
    meta_df.to_csv(RESULTS_DIR / "meta_paths.csv", index=False)

print(f"\nTotal paths found: {len(all_paths_records)}")

# COVID-19 specific: drugs with most path evidence
path_counts = df_paths.groupby("drug_name").size().sort_values(ascending=False)
print("\nDrugs with most COVID-19 paths:")
print(path_counts)

log_event("handoff_completed", {
    "files_written": [
        "results/drug_disease_paths.csv",
        "results/path_narratives.csv",
        "results/meta_paths.csv",
    ],
    "covid_paths_found": len(all_paths_records),
    "status": "ok"
})

print("\n[✓] Path reasoning complete.")
