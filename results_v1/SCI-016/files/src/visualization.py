from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", palette="colorblind")
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 150


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_vdj_annotation(raw: pd.DataFrame, quality: pd.DataFrame, v_usage: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    top_v = v_usage.groupby("v_call")["clone_count"].sum().nlargest(10).index
    heatmap_df = v_usage.loc[v_usage["v_call"].isin(top_v)].pivot(index="sample_id", columns="v_call", values="freq").fillna(0)
    sns.heatmap(heatmap_df, cmap="viridis", ax=axes[0])
    axes[0].set_title("V gene usage")
    axes[0].set_xlabel("V gene")
    axes[0].set_ylabel("Sample")

    productive_lengths = raw.loc[raw["productive"]].copy()
    productive_lengths["cdr3_length"] = productive_lengths["junction_aa"].str.replace("*", "", regex=False).str.len()
    sns.violinplot(data=productive_lengths, x="sample_type", y="cdr3_length", ax=axes[1], palette="colorblind")
    axes[1].set_title("CDR3 length distribution")
    axes[1].set_xlabel("Sample type")
    axes[1].set_ylabel("CDR3 length (aa)")
    axes[1].tick_params(axis="x", rotation=20)

    ratio_df = quality[["sample_id", "productive_ratio"]].copy()
    ratio_df["non_productive_ratio"] = 1 - ratio_df["productive_ratio"]
    ratio_df = ratio_df.set_index("sample_id")
    ratio_df[["productive_ratio", "non_productive_ratio"]].plot(kind="bar", stacked=True, ax=axes[2], color=[sns.color_palette("colorblind")[0], sns.color_palette("colorblind")[3]])
    axes[2].set_title("Productive ratio")
    axes[2].set_xlabel("Sample")
    axes[2].set_ylabel("Read fraction")
    axes[2].legend(title="Status", labels=["Productive", "Non-productive"])
    _save(fig, output_path)


def plot_diversity_metrics(diversity: pd.DataFrame, output_path: Path) -> None:
    metrics = [
        ("shannon_entropy", "Shannon entropy"),
        ("chao1", "Chao1"),
        ("hill_q0", "Hill q=0"),
        ("hill_q1", "Hill q=1"),
        ("hill_q2", "Hill q=2"),
        ("d50", "D50"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    for ax, (col, title) in zip(axes.ravel(), metrics):
        sns.barplot(data=diversity, x="sample_id", y=col, hue="sample_type", dodge=False, ax=ax, palette="colorblind")
        ax.set_title(title)
        ax.set_xlabel("Sample")
        ax.set_ylabel(title)
        ax.tick_params(axis="x", rotation=45)
        if ax is not axes.ravel()[0]:
            ax.legend_.remove()
    axes.ravel()[0].legend(title="Sample type", bbox_to_anchor=(1.05, 1), loc="upper left")
    _save(fig, output_path)


def plot_clonotype_distribution(clonotypes: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    for sample_id, group in clonotypes.groupby("sample_id"):
        ranked = group.sort_values("clone_count", ascending=False)["clone_count"].to_numpy()
        axes[0].plot(np.arange(1, len(ranked) + 1), ranked, label=sample_id)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_title("Rank-abundance curves")
    axes[0].set_xlabel("Clone rank")
    axes[0].set_ylabel("Clone count")
    axes[0].legend(ncol=2, fontsize=8)

    top20 = clonotypes.nlargest(20, "clone_count").copy()
    top20["label"] = top20["v_call"] + "\n" + top20["junction_aa"].str[:8]
    stacked = top20.pivot_table(index="sample_id", columns="label", values="clone_frequency", aggfunc="sum", fill_value=0)
    stacked.plot(kind="bar", stacked=True, ax=axes[1], colormap="viridis")
    axes[1].set_title("Top 20 clonotypes")
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Clone frequency")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=7)
    _save(fig, output_path)


def plot_public_tcr_hla(public_tcrs: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6))
    bubble = (
        public_tcrs.groupby(["antigen", "hla_restriction"], as_index=False)
        .agg(count=("junction_aa", "size"), confidence=("confidence_score", "mean"))
    )
    antigen_order = {name: i for i, name in enumerate(sorted(bubble["antigen"].unique()))}
    hla_order = {name: i for i, name in enumerate(sorted(bubble["hla_restriction"].unique()))}
    scatter = ax.scatter(
        bubble["hla_restriction"].map(hla_order),
        bubble["antigen"].map(antigen_order),
        s=bubble["count"] * 60,
        c=bubble["confidence"],
        cmap="viridis",
        alpha=0.8,
        edgecolor="black",
    )
    ax.set_xticks(list(hla_order.values()), list(hla_order.keys()), rotation=25)
    ax.set_yticks(list(antigen_order.values()), list(antigen_order.keys()))
    ax.set_title("Public TCRs and HLA restriction")
    ax.set_xlabel("HLA type")
    ax.set_ylabel("Antigen")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Mean confidence")
    _save(fig, output_path)


def plot_epitope_models(epitope_metrics_path: Path, output_path: Path) -> None:
    metrics = json.loads(epitope_metrics_path.read_text(encoding="utf-8"))
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 5)
    ax1 = fig.add_subplot(gs[0, :2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    heat_axes = [fig.add_subplot(gs[1, i]) for i in range(5)]

    ax1.plot(metrics["cnn_history"], label="CNN", linewidth=2)
    ax1.plot(metrics["transformer_history"], label="Transformer", linewidth=2)
    ax1.set_title("Training loss curves")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Binary cross-entropy")
    ax1.legend()

    for key, label in [("cnn", "CNN"), ("transformer", "Transformer"), ("ensemble", "Ensemble")]:
        ax2.plot(metrics["roc_curves"][key]["fpr"], metrics["roc_curves"][key]["tpr"], label=f"{label} AUC={metrics[f'{key}_auc']:.3f}")
    ax2.plot([0, 1], [0, 1], linestyle="--", color="gray")
    ax2.set_title("ROC curves")
    ax2.set_xlabel("False positive rate")
    ax2.set_ylabel("True positive rate")
    ax2.legend()

    for ax, matrix, pair in zip(heat_axes, metrics["attention_maps"], metrics["example_pairs"]):
        arr = np.array(matrix)
        sns.heatmap(arr, cmap="viridis", cbar=False, ax=ax)
        ax.set_title(f"{pair['tcr_seq'][:6]} vs {pair['epitope'][:5]}")
        ax.set_xlabel("Token")
        ax.set_ylabel("Token")
    _save(fig, output_path)


def plot_immune_age(immune_age: pd.DataFrame, output_path: Path) -> None:
    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], polar=True)

    sns.scatterplot(
        data=immune_age,
        x="chronological_age",
        y="immune_age_score",
        hue="sample_type",
        style="immunologically_aged",
        s=110,
        palette="colorblind",
        ax=ax1,
    )
    lims = [min(immune_age["chronological_age"].min(), immune_age["immune_age_score"].min()) - 2, max(immune_age["chronological_age"].max(), immune_age["immune_age_score"].max()) + 2]
    ax1.plot(lims, lims, linestyle="--", color="gray")
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_title("Immune age vs chronological age")
    ax1.set_xlabel("Chronological age")
    ax1.set_ylabel("Immune age score")

    radar_features = ["singleton_ratio", "shannon_entropy", "top10_clone_frequency", "public_tcr_count", "mean_cdr3_length"]
    radar = immune_age.groupby("sample_type")[radar_features].mean()
    normalized = (radar - radar.min()) / (radar.max() - radar.min() + 1e-9)
    categories = radar_features
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]
    for sample_type, row in normalized.iterrows():
        values = row.tolist()
        values += values[:1]
        ax2.plot(angles, values, label=sample_type)
        ax2.fill(angles, values, alpha=0.15)
    ax2.set_xticks(angles[:-1], [c.replace("_", " ") for c in categories])
    ax2.set_title("Feature radar chart")
    ax2.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1))
    _save(fig, output_path)


def plot_icb_biomarkers(icb_results: dict, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    metrics = icb_results["metrics"]
    axes[0].plot(metrics["roc_curve"]["fpr"], metrics["roc_curve"]["tpr"], color=sns.color_palette("colorblind")[0], linewidth=2)
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="gray")
    axes[0].set_title(f"ICB response ROC (AUC={metrics['ensemble_auc']:.3f})")
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")

    importance = icb_results["importance"].head(8)
    sns.barplot(data=importance, y="feature", x="importance_mean", ax=axes[1], palette="viridis")
    axes[1].set_title("Feature importance")
    axes[1].set_xlabel("Mean importance")
    axes[1].set_ylabel("Feature")

    embedding = icb_results["embedding"]
    scatter = axes[2].scatter(embedding["x"], embedding["y"], c=embedding["probability"], cmap="viridis", s=110, edgecolor="black")
    for row in embedding.itertuples(index=False):
        axes[2].text(row.x + 0.02, row.y + 0.02, row.sample_id, fontsize=8)
    axes[2].set_title("UMAP-like sample embedding")
    axes[2].set_xlabel("Component 1")
    axes[2].set_ylabel("Component 2")
    cbar = fig.colorbar(scatter, ax=axes[2])
    cbar.set_label("Predicted response probability")
    _save(fig, output_path)
