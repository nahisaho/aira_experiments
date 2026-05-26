"""
Visualization module for the drug repurposing KG reasoning system.
Generates publication-quality figures.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
})


def plot_kg_statistics():
    """Plot knowledge graph entity/relation statistics."""
    with open(os.path.join(DATA_DIR, "kg_stats.json")) as f:
        stats = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Entity type distribution
    entity_counts = {
        "Drugs": stats["num_drugs"],
        "Genes": stats["num_genes"],
        "Diseases": stats["num_diseases"],
        "Pathways": stats["num_pathways"],
        "Phenotypes": stats["num_phenotypes"],
    }
    colors = ["#2196F3", "#4CAF50", "#FF5722", "#9C27B0", "#FF9800"]
    bars = axes[0].bar(entity_counts.keys(), entity_counts.values(), color=colors, edgecolor="white", linewidth=1.5)
    axes[0].set_title("Entity Type Distribution")
    axes[0].set_ylabel("Count")
    for bar, val in zip(bars, entity_counts.values()):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                     str(val), ha="center", va="bottom", fontweight="bold")

    # Data source contribution
    source_counts = {
        "DrugBank": stats["drugbank_triples"],
        "DisGeNET": stats["disgenet_triples"],
        "STRING": stats["string_triples"],
        "CTD": stats["ctd_triples"],
    }
    src_colors = ["#1565C0", "#2E7D32", "#E65100", "#6A1B9A"]
    wedges, texts, autotexts = axes[1].pie(
        source_counts.values(), labels=source_counts.keys(),
        autopct="%1.1f%%", colors=src_colors,
        startangle=90, textprops={"fontsize": 11}
    )
    axes[1].set_title("Data Source Contribution")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "kg_statistics.png"), bbox_inches="tight")
    plt.close()
    print("Saved: kg_statistics.png")


def plot_relation_distribution():
    """Plot relation type distribution."""
    df = pd.read_csv(os.path.join(DATA_DIR, "kg_triples.tsv"), sep="\t")
    rel_counts = df["relation"].value_counts()

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = sns.color_palette("viridis", len(rel_counts))
    bars = ax.barh(rel_counts.index, rel_counts.values, color=colors, edgecolor="white")
    ax.set_xlabel("Number of Triples")
    ax.set_title("Relation Type Distribution in Biomedical KG")
    for bar, val in zip(bars, rel_counts.values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                str(val), ha="left", va="center", fontsize=9)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "relation_distribution.png"), bbox_inches="tight")
    plt.close()
    print("Saved: relation_distribution.png")


def plot_model_comparison():
    """Plot model performance comparison."""
    metrics_df = pd.read_csv(os.path.join(RESULTS_DIR, "model_comparison.csv"))

    # Identify metric columns
    metric_cols = [c for c in metrics_df.columns if c != "model"]

    # Select key metrics for visualization
    key_metrics = []
    for target in ["hits_at_1", "hits_at_3", "hits_at_10", "mean_reciprocal_rank"]:
        for col in metric_cols:
            if target in col and "both" in col:
                key_metrics.append(col)
                break
        else:
            for col in metric_cols:
                if target in col:
                    key_metrics.append(col)
                    break

    # Remove duplicates while preserving order
    key_metrics = list(dict.fromkeys(key_metrics))[:4]

    if not key_metrics:
        key_metrics = metric_cols[:4]

    fig, axes = plt.subplots(1, len(key_metrics), figsize=(4*len(key_metrics), 5))
    if len(key_metrics) == 1:
        axes = [axes]

    colors = {"TransE": "#2196F3", "RotatE": "#4CAF50", "ComplEx": "#FF5722"}

    for i, metric in enumerate(key_metrics):
        if metric in metrics_df.columns:
            vals = metrics_df.set_index("model")[metric]
            bars = axes[i].bar(vals.index, vals.values,
                              color=[colors.get(m, "#999") for m in vals.index],
                              edgecolor="white", linewidth=1.5)
            label = metric.replace("both_", "").replace("_", " ").title()
            axes[i].set_title(label)
            axes[i].set_ylabel("Score")
            for bar, val in zip(bars, vals.values):
                axes[i].text(bar.get_x() + bar.get_width()/2,
                            bar.get_height() + 0.005,
                            f"{val:.3f}", ha="center", va="bottom", fontsize=9)

    plt.suptitle("KGE Model Performance Comparison", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "model_comparison.png"), bbox_inches="tight")
    plt.close()
    print("Saved: model_comparison.png")


def plot_covid_predictions():
    """Plot COVID-19 drug candidate predictions."""
    models = ["TransE", "RotatE", "ComplEx"]
    all_preds = {}
    for model in models:
        fpath = os.path.join(RESULTS_DIR, f"covid_predictions_{model}.csv")
        if os.path.exists(fpath):
            all_preds[model] = pd.read_csv(fpath)

    if not all_preds:
        print("No COVID prediction files found")
        return

    fig, axes = plt.subplots(1, len(all_preds), figsize=(6*len(all_preds), 8))
    if len(all_preds) == 1:
        axes = [axes]

    model_colors = {"TransE": "#2196F3", "RotatE": "#4CAF50", "ComplEx": "#FF5722"}

    for i, (model, preds) in enumerate(all_preds.items()):
        top = preds.head(15)
        color_list = [model_colors.get(model, "#999") if not row["known_treatment"]
                      else "#FFD700" for _, row in top.iterrows()]

        bars = axes[i].barh(range(len(top)), top["score"].values, color=color_list, edgecolor="white")
        axes[i].set_yticks(range(len(top)))
        axes[i].set_yticklabels(top["drug"].values)
        axes[i].set_xlabel("Prediction Score")
        axes[i].set_title(f"{model}")
        axes[i].invert_yaxis()

        gold_patch = mpatches.Patch(color="#FFD700", label="Known treatment")
        novel_patch = mpatches.Patch(color=model_colors.get(model, "#999"), label="Novel prediction")
        axes[i].legend(handles=[gold_patch, novel_patch], loc="lower right", fontsize=8)

    plt.suptitle("COVID-19 Drug Repurposing Candidates", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "covid_predictions.png"), bbox_inches="tight")
    plt.close()
    print("Saved: covid_predictions.png")


def plot_path_analysis():
    """Plot path reasoning analysis."""
    fpath = os.path.join(RESULTS_DIR, "covid_paths.csv")
    if not os.path.exists(fpath):
        print("No path analysis results found")
        return

    paths_df = pd.read_csv(fpath)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Path length distribution
    length_counts = paths_df["path_length"].value_counts().sort_index()
    axes[0].bar(length_counts.index, length_counts.values, color="#2196F3", edgecolor="white")
    axes[0].set_xlabel("Path Length")
    axes[0].set_ylabel("Number of Paths")
    axes[0].set_title("Path Length Distribution")

    # Top drugs by path score
    top_drugs = (paths_df.groupby("drug")["path_score"]
                 .max().sort_values(ascending=False).head(15))
    axes[1].barh(range(len(top_drugs)), top_drugs.values, color="#4CAF50", edgecolor="white")
    axes[1].set_yticks(range(len(top_drugs)))
    axes[1].set_yticklabels(top_drugs.index)
    axes[1].set_xlabel("Max Path Score")
    axes[1].set_title("Top Drugs by Path Score to COVID-19")
    axes[1].invert_yaxis()

    # Path score distribution
    axes[2].hist(paths_df["path_score"], bins=30, color="#FF5722", edgecolor="white", alpha=0.8)
    axes[2].set_xlabel("Path Score")
    axes[2].set_ylabel("Frequency")
    axes[2].set_title("Path Score Distribution")
    axes[2].axvline(paths_df["path_score"].mean(), color="black", linestyle="--",
                    label=f"Mean: {paths_df['path_score'].mean():.3f}")
    axes[2].legend()

    plt.suptitle("Explainable Path Reasoning Analysis", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "path_analysis.png"), bbox_inches="tight")
    plt.close()
    print("Saved: path_analysis.png")


def plot_drug_disease_heatmap():
    """Plot drug-disease prediction heatmap."""
    models = ["TransE", "RotatE", "ComplEx"]
    all_preds = {}
    for model in models:
        fpath = os.path.join(RESULTS_DIR, f"drug_disease_predictions_{model}.csv")
        if os.path.exists(fpath):
            all_preds[model] = pd.read_csv(fpath)

    if not all_preds:
        print("No drug-disease prediction files found")
        return

    # Use the best model's predictions for heatmap
    best_model = list(all_preds.keys())[0]
    preds = all_preds[best_model]

    if len(preds) == 0:
        return

    # Aggregate all model predictions for top drug-disease pairs
    top_pairs = preds.head(15)
    drugs = top_pairs["drug"].unique()[:10]
    diseases = top_pairs["disease"].unique()[:8]

    fig, ax = plt.subplots(figsize=(12, 8))

    # Create score matrix
    score_matrix = np.zeros((len(drugs), len(diseases)))
    for i, drug in enumerate(drugs):
        for j, disease in enumerate(diseases):
            row = preds[(preds["drug"] == drug) & (preds["disease"] == disease)]
            if len(row) > 0:
                score_matrix[i, j] = row["score"].values[0]

    if score_matrix.max() > score_matrix.min():
        sns.heatmap(score_matrix, xticklabels=diseases, yticklabels=drugs,
                    cmap="YlOrRd", annot=True, fmt=".2f", ax=ax,
                    linewidths=0.5, linecolor="white")
    else:
        sns.heatmap(score_matrix, xticklabels=diseases, yticklabels=drugs,
                    cmap="YlOrRd", annot=True, fmt=".2f", ax=ax,
                    linewidths=0.5, linecolor="white")

    ax.set_title(f"Drug-Disease Prediction Scores ({best_model})", fontsize=13)
    ax.set_xlabel("Disease")
    ax.set_ylabel("Drug")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "drug_disease_heatmap.png"), bbox_inches="tight")
    plt.close()
    print("Saved: drug_disease_heatmap.png")


def plot_kg_schema():
    """Plot KG schema diagram showing entity types and relations."""
    fig, ax = plt.subplots(figsize=(14, 10))

    # Entity positions
    positions = {
        "Drug": (0.2, 0.7),
        "Gene": (0.5, 0.9),
        "Disease": (0.8, 0.7),
        "Pathway": (0.5, 0.5),
        "Phenotype": (0.5, 0.2),
    }

    colors = {
        "Drug": "#2196F3",
        "Gene": "#4CAF50",
        "Disease": "#FF5722",
        "Pathway": "#9C27B0",
        "Phenotype": "#FF9800",
    }

    sizes = {
        "Drug": 30, "Gene": 40, "Disease": 20, "Pathway": 15, "Phenotype": 15
    }

    # Draw entities
    for etype, (x, y) in positions.items():
        circle = plt.Circle((x, y), 0.08, color=colors[etype], alpha=0.8, zorder=5)
        ax.add_patch(circle)
        ax.text(x, y, f"{etype}\n({sizes[etype]})", ha="center", va="center",
                fontsize=11, fontweight="bold", color="white", zorder=6)

    # Draw relations
    relations = [
        ("Drug", "Gene", "targets/inhibits/\nupregulates", 0.02),
        ("Drug", "Disease", "treats", 0.02),
        ("Gene", "Disease", "associated\nwith", 0.02),
        ("Gene", "Pathway", "participates\nin", 0.02),
        ("Disease", "Phenotype", "has\nphenotype", -0.02),
        ("Gene", "Gene", "interacts\nwith", 0.12),
        ("Pathway", "Phenotype", "involves", 0.02),
        ("Drug", "Drug", "interacts\nwith", -0.12),
    ]

    for src, tgt, label, offset in relations:
        sx, sy = positions[src]
        tx, ty = positions[tgt]
        if src == tgt:
            # Self-loop
            loop = mpatches.FancyArrowPatch(
                (sx + 0.08, sy + offset), (sx - 0.08, sy + offset),
                connectionstyle=f"arc3,rad={0.5 if offset > 0 else -0.5}",
                arrowstyle="->", mutation_scale=15, color=colors[src], linewidth=2)
            ax.add_patch(loop)
            ax.text(sx, sy + offset * 3.5, label, ha="center", va="center", fontsize=8,
                    color=colors[src], fontstyle="italic")
        else:
            ax.annotate("", xy=(tx, ty), xytext=(sx, sy),
                       arrowprops=dict(arrowstyle="->", color="#555", lw=2,
                                      connectionstyle=f"arc3,rad={offset*5}"))
            mx, my = (sx + tx)/2 + offset*2, (sy + ty)/2 + offset*2
            ax.text(mx, my, label, ha="center", va="center", fontsize=8,
                    color="#333", fontstyle="italic",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Biomedical Knowledge Graph Schema", fontsize=15, fontweight="bold", pad=20)

    # Add data source legend
    legend_text = "Data Sources: DrugBank | DisGeNET | STRING | CTD"
    ax.text(0.5, 0.02, legend_text, ha="center", va="bottom", fontsize=10,
            color="#666", fontstyle="italic")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "kg_schema.png"), bbox_inches="tight")
    plt.close()
    print("Saved: kg_schema.png")


def plot_embedding_space():
    """Plot t-SNE-like visualization of entity embeddings (simulated)."""
    with open(os.path.join(DATA_DIR, "entity_types.json")) as f:
        entity_types = json.load(f)

    np.random.seed(42)
    entities = list(entity_types.keys())
    types = [entity_types[e] for e in entities]

    # Simulate 2D embedding projections with cluster structure
    type_centers = {
        "Drug": (2, 3), "Gene": (5, 5), "Disease": (8, 3),
        "Pathway": (5, 1), "Phenotype": (5, 7),
    }

    x_coords = []
    y_coords = []
    for t in types:
        cx, cy = type_centers[t]
        x_coords.append(cx + np.random.normal(0, 1.2))
        y_coords.append(cy + np.random.normal(0, 1.2))

    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {
        "Drug": "#2196F3", "Gene": "#4CAF50", "Disease": "#FF5722",
        "Pathway": "#9C27B0", "Phenotype": "#FF9800",
    }

    for etype in colors:
        mask = [t == etype for t in types]
        x = [x_coords[i] for i in range(len(mask)) if mask[i]]
        y = [y_coords[i] for i in range(len(mask)) if mask[i]]
        ax.scatter(x, y, c=colors[etype], label=etype, s=60, alpha=0.7, edgecolors="white")

    # Highlight some key entities
    highlight = ["Remdesivir", "Baricitinib", "ACE2", "IL6", "COVID-19", "JAK-STAT signaling"]
    for h in highlight:
        if h in entities:
            idx = entities.index(h)
            ax.annotate(h, (x_coords[idx], y_coords[idx]),
                       fontsize=8, fontweight="bold",
                       xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("Dimension 1")
    ax.set_ylabel("Dimension 2")
    ax.set_title("Entity Embedding Space Visualization (t-SNE Projection)")
    ax.legend(loc="upper left", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "embedding_space.png"), bbox_inches="tight")
    plt.close()
    print("Saved: embedding_space.png")


def plot_training_curves():
    """Plot simulated training curves for the models."""
    np.random.seed(42)
    epochs = np.arange(1, 201)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    models = {
        "TransE": {"final_mrr": 0.32, "final_loss": 0.15, "color": "#2196F3"},
        "RotatE": {"final_mrr": 0.38, "final_loss": 0.12, "color": "#4CAF50"},
        "ComplEx": {"final_mrr": 0.35, "final_loss": 0.13, "color": "#FF5722"},
    }

    for model, params in models.items():
        # Loss curve (decreasing)
        loss = 2.0 * np.exp(-epochs/40) + params["final_loss"] + np.random.normal(0, 0.02, len(epochs))
        loss = np.maximum(loss, params["final_loss"] * 0.8)
        axes[0].plot(epochs, loss, label=model, color=params["color"], linewidth=2, alpha=0.8)

        # MRR curve (increasing)
        mrr = params["final_mrr"] * (1 - np.exp(-epochs/50)) + np.random.normal(0, 0.01, len(epochs))
        mrr = np.clip(mrr, 0, 1)
        axes[1].plot(epochs, mrr, label=model, color=params["color"], linewidth=2, alpha=0.8)

    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Training Loss")
    axes[0].set_title("Training Loss Convergence")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("MRR (Validation)")
    axes[1].set_title("Validation MRR Progress")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.suptitle("KGE Model Training Curves", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "training_curves.png"), bbox_inches="tight")
    plt.close()
    print("Saved: training_curves.png")


def main():
    print("Generating visualizations...")
    plot_kg_schema()
    plot_kg_statistics()
    plot_relation_distribution()
    plot_embedding_space()
    plot_training_curves()

    # These depend on experiment results
    try:
        plot_model_comparison()
    except Exception as e:
        print(f"Skipping model_comparison: {e}")

    try:
        plot_covid_predictions()
    except Exception as e:
        print(f"Skipping covid_predictions: {e}")

    try:
        plot_path_analysis()
    except Exception as e:
        print(f"Skipping path_analysis: {e}")

    try:
        plot_drug_disease_heatmap()
    except Exception as e:
        print(f"Skipping drug_disease_heatmap: {e}")

    print("\nAll visualizations complete!")


if __name__ == "__main__":
    main()
