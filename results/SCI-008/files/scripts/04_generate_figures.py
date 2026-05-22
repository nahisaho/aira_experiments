#!/usr/bin/env python3
"""
Step 4: Generate all figures for the report and paper.
"""

import json
import os
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
LOG_DIR = os.path.join(BASE_DIR, "logs")
os.makedirs(FIGURES_DIR, exist_ok=True)


def fig1_kg_schema():
    """Knowledge Graph schema diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    positions = {
        "Drug": (0.2, 0.7),
        "Gene": (0.5, 0.5),
        "Disease": (0.8, 0.7),
        "Pathway": (0.5, 0.2),
        "Phenotype": (0.8, 0.3),
    }

    colors = {
        "Drug": "#2196F3",
        "Gene": "#4CAF50",
        "Disease": "#F44336",
        "Pathway": "#FF9800",
        "Phenotype": "#9C27B0",
    }

    with open(os.path.join(DATA_DIR, "kg_stats.json")) as f:
        stats = json.load(f)

    for entity, (x, y) in positions.items():
        count = stats["entity_types"].get(entity, 0)
        circle = plt.Circle((x, y), 0.08, color=colors[entity], alpha=0.8)
        ax.add_patch(circle)
        ax.text(x, y + 0.01, entity, ha="center", va="center", fontsize=12, fontweight="bold", color="white")
        ax.text(x, y - 0.03, f"n={count}", ha="center", va="center", fontsize=9, color="white")

    edges = [
        ("Drug", "Gene", "targets"),
        ("Gene", "Disease", "associated_with"),
        ("Gene", "Gene", "interacts_with"),
        ("Gene", "Pathway", "participates_in"),
        ("Disease", "Phenotype", "has_phenotype"),
        ("Drug", "Disease", "treats"),
    ]

    for src, tgt, rel in edges:
        sx, sy = positions[src]
        tx, ty = positions[tgt]
        if src == tgt:
            ax.annotate("", xy=(sx + 0.06, sy + 0.06), xytext=(sx - 0.06, sy + 0.06),
                        arrowprops=dict(arrowstyle="->", color="gray", lw=1.5,
                                       connectionstyle="arc3,rad=0.5"))
            ax.text(sx, sy + 0.14, rel, ha="center", fontsize=8, color="gray", style="italic")
        else:
            mx, my = (sx + tx) / 2, (sy + ty) / 2
            ax.annotate("", xy=(tx, ty), xytext=(sx, sy),
                        arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
            ax.text(mx, my + 0.03, rel, ha="center", fontsize=8, color="gray", style="italic")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Biomedical Knowledge Graph Schema", fontsize=14, fontweight="bold", pad=20)

    total = stats["total_triples"]
    ent = stats["total_entities"]
    ax.text(0.5, 0.95, f"Total: {ent} entities, {total} triples", ha="center",
            fontsize=10, transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_kg_schema.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig1_kg_schema.png")


def fig2_entity_distribution():
    """Entity and relation distribution."""
    with open(os.path.join(DATA_DIR, "kg_stats.json")) as f:
        stats = json.load(f)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Entity type distribution
    types = list(stats["entity_types"].keys())
    counts = list(stats["entity_types"].values())
    colors = ["#2196F3", "#4CAF50", "#F44336", "#FF9800", "#9C27B0"]

    bars1 = ax1.bar(types, counts, color=colors, edgecolor="white", linewidth=1.5)
    ax1.set_title("Entity Type Distribution", fontsize=13, fontweight="bold")
    ax1.set_ylabel("Count", fontsize=11)
    ax1.set_xlabel("Entity Type", fontsize=11)
    for bar, count in zip(bars1, counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 str(count), ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Relation distribution
    rels = list(stats["relation_types"].keys())
    rel_counts = list(stats["relation_types"].values())
    colors2 = sns.color_palette("viridis", len(rels))

    bars2 = ax2.barh(rels, rel_counts, color=colors2, edgecolor="white", linewidth=1.5)
    ax2.set_title("Relation Type Distribution", fontsize=13, fontweight="bold")
    ax2.set_xlabel("Count", fontsize=11)
    for bar, count in zip(bars2, rel_counts):
        ax2.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                 str(count), ha="left", va="center", fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig2_entity_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig2_entity_distribution.png")


def fig3_model_comparison():
    """Model comparison metrics."""
    comp_path = os.path.join(RESULTS_DIR, "model_comparison.csv")
    if not os.path.exists(comp_path):
        print("  Skipping fig3: model_comparison.csv not found")
        return

    df = pd.read_csv(comp_path)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Hits@K comparison
    metrics_hits = ["hits_at_1", "hits_at_3", "hits_at_10"]
    labels_hits = ["Hits@1", "Hits@3", "Hits@10"]
    x = np.arange(len(labels_hits))
    width = 0.25
    colors = ["#2196F3", "#4CAF50", "#F44336"]

    for i, (_, row) in enumerate(df.iterrows()):
        vals = [row.get(m, 0) or 0 for m in metrics_hits]
        axes[0].bar(x + i * width, vals, width, label=row["model"], color=colors[i], alpha=0.85)

    axes[0].set_xticks(x + width)
    axes[0].set_xticklabels(labels_hits)
    axes[0].set_ylabel("Score")
    axes[0].set_title("Hits@K Comparison", fontweight="bold")
    axes[0].legend()
    axes[0].set_ylim(0, 1.0)

    # MRR comparison
    mrr_vals = df["mean_reciprocal_rank"].fillna(0).values
    axes[1].bar(df["model"], mrr_vals, color=colors, alpha=0.85)
    axes[1].set_ylabel("MRR")
    axes[1].set_title("Mean Reciprocal Rank", fontweight="bold")
    axes[1].set_ylim(0, max(mrr_vals) * 1.3 if max(mrr_vals) > 0 else 1.0)
    for i, v in enumerate(mrr_vals):
        axes[1].text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=10)

    # Training time
    times = df["training_time_sec"].values
    axes[2].bar(df["model"], times, color=colors, alpha=0.85)
    axes[2].set_ylabel("Time (seconds)")
    axes[2].set_title("Training Time", fontweight="bold")
    for i, v in enumerate(times):
        axes[2].text(i, v + 1, f"{v:.0f}s", ha="center", fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig3_model_comparison.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig3_model_comparison.png")


def fig4_covid_predictions():
    """COVID-19 drug repurposing predictions."""
    pred_path = os.path.join(RESULTS_DIR, "covid19_predictions.csv")
    if not os.path.exists(pred_path):
        print("  Skipping fig4: covid19_predictions.csv not found")
        return

    df = pd.read_csv(pred_path)
    top20 = df.head(20)

    fig, ax = plt.subplots(figsize=(12, 8))

    colors = ["#F44336" if k else "#2196F3" for k in top20["known"]]
    bars = ax.barh(range(len(top20)), top20["score"], color=colors, alpha=0.85, edgecolor="white")

    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20["drug_name"], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Prediction Score", fontsize=12)
    ax.set_title("Top 20 Predicted Drugs for COVID-19", fontsize=14, fontweight="bold")

    known_patch = mpatches.Patch(color="#F44336", label="Known COVID-19 Drug")
    novel_patch = mpatches.Patch(color="#2196F3", label="Novel Prediction")
    ax.legend(handles=[known_patch, novel_patch], loc="lower right", fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig4_covid_predictions.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig4_covid_predictions.png")


def fig5_path_explanation():
    """Explainable path reasoning visualization."""
    expl_path = os.path.join(RESULTS_DIR, "covid19_path_explanations.json")
    if not os.path.exists(expl_path):
        print("  Skipping fig5: covid19_path_explanations.json not found")
        return

    with open(expl_path) as f:
        explanations = json.load(f)

    # Find first drug with paths
    drug_name = None
    paths = []
    for dn, p in explanations.items():
        if len(p) > 0:
            drug_name = dn
            paths = p[:5]
            break

    if not drug_name or not paths:
        print("  Skipping fig5: no paths found")
        return

    fig, ax = plt.subplots(figsize=(14, 8))

    y_offset = 0.9
    ax.text(0.5, 0.98, f"Explanatory Paths: {drug_name} → COVID-19",
            ha="center", va="top", fontsize=14, fontweight="bold", transform=ax.transAxes)

    type_colors = {
        "targets": "#2196F3",
        "associated_with": "#F44336",
        "interacts_with": "#4CAF50",
        "participates_in": "#FF9800",
        "has_phenotype": "#9C27B0",
        "treats": "#795548",
    }

    for pi, path in enumerate(paths):
        y = y_offset - pi * 0.18
        x_start = 0.05
        x_step = 0.85 / max(len(path), 1)

        for si, step in enumerate(path):
            x = x_start + si * x_step
            # Node
            ax.text(x, y, step["from_name"], fontsize=8, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#E3F2FD", edgecolor="#1565C0"))
            # Arrow + relation
            rel_color = type_colors.get(step["relation"], "gray")
            ax.annotate("", xy=(x + x_step * 0.6, y), xytext=(x + x_step * 0.15, y),
                        arrowprops=dict(arrowstyle="->", color=rel_color, lw=2))
            ax.text(x + x_step * 0.38, y + 0.04, step["relation"], fontsize=7,
                    ha="center", color=rel_color, style="italic")

            # Last node
            if si == len(path) - 1:
                ax.text(x + x_step * 0.75, y, step["to_name"], fontsize=8, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFEBEE", edgecolor="#C62828"))

        ax.text(0.01, y, f"P{pi+1}", fontsize=9, fontweight="bold", va="center", color="gray")

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig5_path_explanation.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig5_path_explanation.png")


def fig6_kg_subgraph():
    """COVID-19-centric KG subgraph visualization."""
    df = pd.read_csv(os.path.join(DATA_DIR, "triples.tsv"), sep="\t")
    with open(os.path.join(DATA_DIR, "entities.json")) as f:
        entities = json.load(f)

    # Build subgraph around COVID-19
    covid_id = "DOID:0080600"
    G = nx.DiGraph()

    # Nodes connected to COVID-19
    covid_neighbors = set()
    for _, row in df.iterrows():
        if row["head"] == covid_id or row["tail"] == covid_id:
            covid_neighbors.add(row["head"])
            covid_neighbors.add(row["tail"])
            G.add_edge(row["head"], row["tail"], relation=row["relation"])

    # 2-hop neighbors (limited)
    for _, row in df.iterrows():
        if row["head"] in covid_neighbors or row["tail"] in covid_neighbors:
            if row["head"] in covid_neighbors and row["tail"] in covid_neighbors:
                G.add_edge(row["head"], row["tail"], relation=row["relation"])

    # Add select drug connections
    for _, row in df.iterrows():
        if row["relation"] == "treats" and row["tail"] == covid_id:
            G.add_edge(row["head"], row["tail"], relation=row["relation"])
            # Add targets of these drugs
            drug_targets = df[(df["head"] == row["head"]) & (df["relation"] == "targets")]
            for _, tr in drug_targets.iterrows():
                if tr["tail"] in covid_neighbors:
                    G.add_edge(tr["head"], tr["tail"], relation=tr["relation"])

    fig, ax = plt.subplots(figsize=(16, 12))

    type_colors = {}
    for node in G.nodes():
        info = entities.get(node, {})
        t = info.get("type", "Unknown")
        if t == "Drug":
            type_colors[node] = "#2196F3"
        elif t == "Gene":
            type_colors[node] = "#4CAF50"
        elif t == "Disease":
            type_colors[node] = "#F44336"
        elif t == "Pathway":
            type_colors[node] = "#FF9800"
        elif t == "Phenotype":
            type_colors[node] = "#9C27B0"
        else:
            type_colors[node] = "#9E9E9E"

    node_colors = [type_colors.get(n, "#9E9E9E") for n in G.nodes()]
    node_sizes = [800 if n == covid_id else 300 for n in G.nodes()]
    labels = {n: entities.get(n, {}).get("name", n)[:15] for n in G.nodes()}

    pos = nx.spring_layout(G, seed=42, k=2)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color="#BDBDBD", alpha=0.4, arrows=True,
                           arrowsize=10, ax=ax)

    legend_elements = [
        mpatches.Patch(facecolor="#2196F3", label="Drug"),
        mpatches.Patch(facecolor="#4CAF50", label="Gene"),
        mpatches.Patch(facecolor="#F44336", label="Disease"),
        mpatches.Patch(facecolor="#FF9800", label="Pathway"),
        mpatches.Patch(facecolor="#9C27B0", label="Phenotype"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=10)
    ax.set_title("COVID-19-Centric Knowledge Graph Subgraph", fontsize=14, fontweight="bold")
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig6_kg_subgraph.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig6_kg_subgraph.png")


def fig7_heatmap_drug_disease():
    """Drug-Disease prediction heatmap for COVID-related diseases."""
    pred_path = os.path.join(RESULTS_DIR, "all_drug_disease_predictions.csv")
    if not os.path.exists(pred_path):
        print("  Skipping fig7: all_drug_disease_predictions.csv not found")
        return

    df = pd.read_csv(pred_path)

    # Focus on COVID-related diseases and top drugs
    target_diseases = [
        "COVID-19", "Cytokine Storm", "Acute Respiratory Distress",
        "Pneumonia", "Thrombotic Disorder",
    ]
    top_drugs = df[df["disease_name"].isin(target_diseases)].groupby("drug_name")["score"].mean()
    top_drugs = top_drugs.nlargest(15).index.tolist()

    subset = df[df["drug_name"].isin(top_drugs) & df["disease_name"].isin(target_diseases)]
    if len(subset) == 0:
        print("  Skipping fig7: no matching predictions")
        return

    pivot = subset.pivot_table(index="drug_name", columns="disease_name", values="score", aggfunc="first")

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlOrRd", linewidths=0.5, ax=ax,
                cbar_kws={"label": "Prediction Score"})
    ax.set_title("Drug-Disease Prediction Scores (COVID-Related)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Drug")
    ax.set_xlabel("Disease")
    plt.xticks(rotation=30, ha="right")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig7_heatmap_drug_disease.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig7_heatmap_drug_disease.png")


def fig8_degree_distribution():
    """Degree distribution of the KG."""
    df = pd.read_csv(os.path.join(DATA_DIR, "triples.tsv"), sep="\t")
    G = nx.DiGraph()
    for _, row in df.iterrows():
        G.add_edge(row["head"], row["tail"])

    degrees = [d for _, d in G.degree()]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.hist(degrees, bins=30, color="#2196F3", alpha=0.8, edgecolor="white")
    ax1.set_xlabel("Degree", fontsize=11)
    ax1.set_ylabel("Frequency", fontsize=11)
    ax1.set_title("Degree Distribution", fontsize=13, fontweight="bold")
    ax1.axvline(np.mean(degrees), color="red", linestyle="--", label=f"Mean={np.mean(degrees):.1f}")
    ax1.legend()

    # Log-log
    from collections import Counter
    deg_count = Counter(degrees)
    degs = sorted(deg_count.keys())
    counts = [deg_count[d] for d in degs]
    ax2.scatter(degs, counts, color="#4CAF50", alpha=0.7, s=30)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Degree (log)", fontsize=11)
    ax2.set_ylabel("Frequency (log)", fontsize=11)
    ax2.set_title("Degree Distribution (Log-Log)", fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig8_degree_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig8_degree_distribution.png")


def main():
    print("=== Generating Figures ===")
    fig1_kg_schema()
    fig2_entity_distribution()
    fig3_model_comparison()
    fig4_covid_predictions()
    fig5_path_explanation()
    fig6_kg_subgraph()
    fig7_heatmap_drug_disease()
    fig8_degree_distribution()

    # Log
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": "visualization",
        "event_type": "figure_generation",
        "actor": "co-scientist",
        "skill_or_tool": "04_generate_figures.py",
        "files_written": [f"figures/{f}" for f in os.listdir(FIGURES_DIR) if f.endswith(".png")],
        "status": "ok",
    }
    with open(os.path.join(LOG_DIR, "process-log.jsonl"), "a") as f:
        f.write(json.dumps(log_entry) + "\n")

    print("=== Figure generation complete ===")


if __name__ == "__main__":
    main()
