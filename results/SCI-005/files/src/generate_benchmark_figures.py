from __future__ import annotations

import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
SEED = 42


def setup_environment() -> None:
    random.seed(SEED)
    np.random.seed(SEED)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk", palette="colorblind")
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "legend.frameon": True,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    png_path = FIGURES_DIR / f"{stem}.png"
    svg_path = FIGURES_DIR / f"{stem}.svg"
    fig.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


SV_PERFORMANCE = {
    "DeepSV-LR": {
        "DEL": (0.96, 0.94, 0.95),
        "INS": (0.94, 0.92, 0.93),
        "DUP": (0.91, 0.88, 0.895),
        "INV": (0.89, 0.85, 0.87),
        "BND/TRA": (0.87, 0.82, 0.845),
    },
    "Sniffles2": {
        "DEL": (0.94, 0.92, 0.93),
        "INS": (0.92, 0.90, 0.91),
        "DUP": (0.88, 0.84, 0.86),
        "INV": (0.85, 0.80, 0.825),
        "BND/TRA": (0.83, 0.78, 0.805),
    },
    "CuteSV": {
        "DEL": (0.93, 0.91, 0.92),
        "INS": (0.91, 0.89, 0.90),
        "DUP": (0.87, 0.83, 0.85),
        "INV": (0.84, 0.79, 0.815),
        "BND/TRA": (0.82, 0.76, 0.79),
    },
    "SVIM": {
        "DEL": (0.91, 0.89, 0.90),
        "INS": (0.89, 0.87, 0.88),
        "DUP": (0.85, 0.81, 0.83),
        "INV": (0.82, 0.77, 0.795),
        "BND/TRA": (0.80, 0.74, 0.77),
    },
    "pbsv": {
        "DEL": (0.92, 0.90, 0.91),
        "INS": (0.90, 0.88, 0.89),
        "DUP": (0.86, 0.82, 0.84),
        "INV": (0.83, 0.78, 0.805),
        "BND/TRA": (0.81, 0.75, 0.78),
    },
    "DELLY2": {
        "DEL": (0.88, 0.85, 0.865),
        "INS": (0.82, 0.78, 0.80),
        "DUP": (0.80, 0.75, 0.775),
        "INV": (0.78, 0.72, 0.75),
        "BND/TRA": (0.75, 0.70, 0.725),
    },
}


def build_performance_dataframe() -> pd.DataFrame:
    rows = []
    for tool, sv_map in SV_PERFORMANCE.items():
        for sv_type, (precision, recall, f1) in sv_map.items():
            rows.append(
                {
                    "Tool": tool,
                    "SV Type": sv_type,
                    "Precision": precision,
                    "Recall": recall,
                    "F1": f1,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "sv_performance_metrics.csv", index=False)
    return df


def figure_sv_performance(df: pd.DataFrame) -> None:
    metric_palette = sns.color_palette("colorblind", 3)
    melted = df.melt(
        id_vars=["Tool", "SV Type"],
        value_vars=["Precision", "Recall", "F1"],
        var_name="Metric",
        value_name="Score",
    )
    sv_types = ["DEL", "INS", "DUP", "INV", "BND/TRA"]
    fig, axes = plt.subplots(3, 2, figsize=(18, 14), sharey=True)
    axes = axes.flatten()

    for idx, sv_type in enumerate(sv_types):
        ax = axes[idx]
        subset = melted[melted["SV Type"] == sv_type]
        sns.barplot(
            data=subset,
            x="Tool",
            y="Score",
            hue="Metric",
            palette=metric_palette,
            edgecolor="black",
            linewidth=0.4,
            ax=ax,
        )
        ax.set_title(f"{sv_type} detection")
        ax.set_xlabel("")
        ax.set_ylabel("Score")
        ax.set_ylim(0.65, 1.0)
        ax.tick_params(axis="x", rotation=35)
        for label in ax.get_xticklabels():
            if label.get_text() == "DeepSV-LR":
                label.set_fontweight("bold")
        if idx == 0:
            ax.legend(title="Metric", loc="upper right")
        else:
            ax.get_legend().remove()

    axes[-1].axis("off")
    fig.suptitle("SV Detection Performance Comparison", fontsize=22, fontweight="bold", y=1.02)
    save_figure(fig, "sv_performance_comparison")


SENSITIVITY_PARAMS = {
    "DeepSV-LR": {"peak": 0.968, "small_drop": 0.11, "large_drop": 0.045},
    "Sniffles2": {"peak": 0.948, "small_drop": 0.16, "large_drop": 0.085},
    "CuteSV": {"peak": 0.94, "small_drop": 0.19, "large_drop": 0.10},
}


def sensitivity_curve(log_sizes: np.ndarray, peak: float, small_drop: float, large_drop: float) -> np.ndarray:
    small_penalty = small_drop / (1 + np.exp((log_sizes - 2.4) * 4.2))
    large_penalty = large_drop / (1 + np.exp((6.1 - log_sizes) * 3.1))
    mid_sculpt = 0.012 * ((log_sizes - 4.4) ** 2) / 5.5
    curve = peak - small_penalty - large_penalty - mid_sculpt
    return np.clip(curve, 0.45, 0.985)


def figure_size_sensitivity() -> None:
    sizes = np.logspace(np.log10(50), np.log10(1e7), 240)
    log_sizes = np.log10(sizes)
    curve_df = pd.DataFrame({"SV Size (bp)": sizes})

    fig, ax = plt.subplots(figsize=(11, 7))
    palette = sns.color_palette("colorblind", 3)

    for color, (tool, params) in zip(palette, SENSITIVITY_PARAMS.items()):
        sensitivity = sensitivity_curve(log_sizes, **params)
        curve_df[tool] = sensitivity
        ax.plot(sizes, sensitivity, linewidth=3, color=color, label=tool)

    curve_df.to_csv(RESULTS_DIR / "sv_size_sensitivity.csv", index=False)
    ax.axvspan(50, 300, color="#d9d9d9", alpha=0.18)
    ax.axvspan(1e6, 1e7, color="#d9d9d9", alpha=0.18)
    ax.set_xscale("log")
    ax.set_ylim(0.45, 1.0)
    ax.set_xlabel("SV size (bp)")
    ax.set_ylabel("Detection sensitivity")
    ax.set_title("Sensitivity Across Structural Variant Sizes")
    ax.legend(title="Tool", loc="lower right")
    save_figure(fig, "sv_size_sensitivity")


PR_TARGETS = {
    "DEL": {"target_auc": 0.97, "gamma": 2.0},
    "INS": {"target_auc": 0.95, "gamma": 1.7},
    "DUP": {"target_auc": 0.93, "gamma": 1.4},
    "INV": {"target_auc": 0.91, "gamma": 1.2},
}


def make_pr_curve(target_auc: float, gamma: float) -> tuple[np.ndarray, np.ndarray, float]:
    recall = np.linspace(0, 1, 400)
    coefficient = (1 - target_auc) * (gamma + 1)
    precision = 1 - coefficient * np.power(recall, gamma)
    precision = np.clip(precision, 0.55, 1.0)
    auc = np.trapezoid(precision, recall)
    return recall, precision, auc


def figure_precision_recall() -> None:
    fig, ax = plt.subplots(figsize=(10, 7))
    palette = sns.color_palette("viridis", len(PR_TARGETS))
    output_rows = []

    for color, (sv_type, params) in zip(palette, PR_TARGETS.items()):
        recall, precision, auc = make_pr_curve(**params)
        output_rows.extend(
            {
                "SV Type": sv_type,
                "Recall": r,
                "Precision": p,
                "AUC": auc,
            }
            for r, p in zip(recall, precision)
        )
        ax.plot(recall, precision, linewidth=3, color=color, label=f"{sv_type} (AUC={auc:.2f})")

    pd.DataFrame(output_rows).to_csv(RESULTS_DIR / "precision_recall_curves.csv", index=False)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.02)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("DeepSV-LR Precision-Recall Curves")
    ax.legend(title="SV type", loc="lower left")
    save_figure(fig, "precision_recall_curves")


REPEAT_REGION_DATA = pd.DataFrame(
    [
        [0.95, 0.90, 0.88, 0.84, 0.79, 0.74],
        [0.93, 0.85, 0.81, 0.76, 0.70, 0.64],
        [0.91, 0.83, 0.79, 0.73, 0.67, 0.61],
        [0.89, 0.80, 0.75, 0.70, 0.64, 0.58],
        [0.90, 0.81, 0.77, 0.71, 0.65, 0.59],
        [0.86, 0.74, 0.69, 0.63, 0.57, 0.52],
    ],
    index=["DeepSV-LR", "Sniffles2", "CuteSV", "SVIM", "pbsv", "DELLY2"],
    columns=["Non-repeat", "Simple Repeat", "SINE/LINE", "Segmental Dup", "Telomere", "Centromere"],
)


def figure_repeat_region() -> None:
    REPEAT_REGION_DATA.to_csv(RESULTS_DIR / "repeat_region_performance.csv")
    fig, ax = plt.subplots(figsize=(11, 7))
    sns.heatmap(
        REPEAT_REGION_DATA,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        linewidths=0.5,
        cbar_kws={"label": "F1 score"},
        ax=ax,
    )
    ax.set_title("Repeat Region SV Detection Performance")
    ax.set_xlabel("Genomic region")
    ax.set_ylabel("Tool")
    save_figure(fig, "repeat_region_performance")


COMPLEX_SV_DATA = pd.DataFrame(
    {
        "Complex SV": [
            "Chromothripsis",
            "ecDNA",
            "Nested SV",
            "Multi-breakpoint",
            "Reciprocal translocation",
        ],
        "DeepSV-LR": [0.86, 0.79, 0.83, 0.81, 0.78],
        "Sniffles2": [0.80, 0.72, 0.76, 0.74, 0.70],
        "CuteSV": [0.77, 0.68, 0.72, 0.70, 0.66],
    }
)


def figure_complex_sv() -> None:
    COMPLEX_SV_DATA.to_csv(RESULTS_DIR / "complex_sv_detection.csv", index=False)
    melted = COMPLEX_SV_DATA.melt(id_vars="Complex SV", var_name="Tool", value_name="Detection rate")
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.barplot(
        data=melted,
        x="Complex SV",
        y="Detection rate",
        hue="Tool",
        palette=sns.color_palette("colorblind", 3),
        edgecolor="black",
        linewidth=0.4,
        ax=ax,
    )
    ax.set_ylim(0.55, 0.9)
    ax.set_xlabel("")
    ax.set_ylabel("Detection rate")
    ax.set_title("Complex Structural Variant Detection")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(title="Tool", loc="upper right")
    save_figure(fig, "complex_sv_detection")


HYBRID_DATA = pd.DataFrame(
    {
        "Analysis mode": ["Long-read only", "Hybrid (LR+SR)"],
        "Precision": [0.939, 0.952],
        "Recall": [0.918, 0.938],
        "F1": [0.928, 0.945],
        "Breakpoint accuracy (bp)": [28, 17],
    }
)


def figure_hybrid_improvement() -> None:
    HYBRID_DATA.to_csv(RESULTS_DIR / "hybrid_improvement.csv", index=False)
    score_df = HYBRID_DATA.melt(
        id_vars="Analysis mode",
        value_vars=["Precision", "Recall", "F1"],
        var_name="Metric",
        value_name="Score",
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [3.4, 1.6]})
    sns.barplot(
        data=score_df,
        x="Metric",
        y="Score",
        hue="Analysis mode",
        palette=sns.color_palette("colorblind", 2),
        edgecolor="black",
        linewidth=0.4,
        ax=axes[0],
    )
    axes[0].set_ylim(0.85, 0.98)
    axes[0].set_title("DeepSV-LR core metrics")
    axes[0].set_ylabel("Score")
    axes[0].set_xlabel("")
    axes[0].legend(title="Analysis mode", loc="upper left")

    sns.barplot(
        data=HYBRID_DATA,
        x="Analysis mode",
        y="Breakpoint accuracy (bp)",
        hue="Analysis mode",
        palette=sns.color_palette("colorblind", 2),
        dodge=False,
        legend=False,
        edgecolor="black",
        linewidth=0.4,
        ax=axes[1],
    )
    axes[1].set_title("Breakpoint accuracy")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("bp (lower is better)")
    axes[1].tick_params(axis="x", rotation=20)

    fig.suptitle("Hybrid vs Long-read Only Analysis", fontsize=20, fontweight="bold", y=1.02)
    save_figure(fig, "hybrid_improvement")


def write_manifest() -> None:
    manifest = {
        "seed": SEED,
        "figures": [
            "figures/sv_performance_comparison.png",
            "figures/sv_size_sensitivity.png",
            "figures/precision_recall_curves.png",
            "figures/repeat_region_performance.png",
            "figures/complex_sv_detection.png",
            "figures/hybrid_improvement.png",
        ],
        "vector_figures": [
            "figures/sv_performance_comparison.svg",
            "figures/sv_size_sensitivity.svg",
            "figures/precision_recall_curves.svg",
            "figures/repeat_region_performance.svg",
            "figures/complex_sv_detection.svg",
            "figures/hybrid_improvement.svg",
        ],
        "results": [
            "results/sv_performance_metrics.csv",
            "results/sv_size_sensitivity.csv",
            "results/precision_recall_curves.csv",
            "results/repeat_region_performance.csv",
            "results/complex_sv_detection.csv",
            "results/hybrid_improvement.csv",
        ],
        "description": "Simulated benchmark outputs for DeepSV-LR structural variant evaluation figures.",
    }
    (RESULTS_DIR / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    setup_environment()
    performance_df = build_performance_dataframe()
    figure_sv_performance(performance_df)
    figure_size_sensitivity()
    figure_precision_recall()
    figure_repeat_region()
    figure_complex_sv()
    figure_hybrid_improvement()
    write_manifest()
    print("Generated benchmark figures and result tables.")


if __name__ == "__main__":
    main()
