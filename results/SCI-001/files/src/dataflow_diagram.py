"""
Data flow diagram for CRISPR-Cas9 off-target prediction pipeline.
Uses matplotlib patches to draw a publication-quality flowchart.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from pathlib import Path


def draw_box(ax, x, y, w, h, text, color, fontsize=9, text_color="white"):
    rect = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04",
        facecolor=color, edgecolor="white", linewidth=1.5, zorder=3,
    )
    ax.add_patch(rect)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, color=text_color,
            fontweight="bold", zorder=4, wrap=True,
            multialignment="center")


def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="-|>", color="#555555",
            lw=1.5, mutation_scale=15,
        ),
        zorder=2,
    )


def create_dataflow_diagram(save_path: str = "figures/dataflow_diagram.png"):
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)

    COLORS = {
        "input":    "#0072B2",
        "process":  "#009E73",
        "feature":  "#E69F00",
        "model":    "#CC79A7",
        "output":   "#D55E00",
        "eval":     "#56B4E9",
    }

    fig, ax = plt.subplots(figsize=(18, 10))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_facecolor("#F8F8F8")
    fig.patch.set_facecolor("#F8F8F8")

    ax.set_title(
        "CRISPR-Cas9 Off-Target Prediction Pipeline — Data Flow",
        fontsize=15, fontweight="bold", pad=12,
    )

    # ── Row 1: Raw Inputs ─────────────────────────────────────────────────────
    draw_box(ax, 3.0, 9.0, 3.2, 0.8, "GUIDE-seq / CIRCLE-seq\nRaw Data", COLORS["input"])
    draw_box(ax, 9.0, 9.0, 3.2, 0.8, "Reference Genome\n(hg38)", COLORS["input"])
    draw_box(ax, 15.0, 9.0, 3.2, 0.8, "Epigenomics\n(ATAC-seq, WGBS)", COLORS["input"])

    # ── Row 2: Preprocessing ─────────────────────────────────────────────────
    draw_box(ax, 3.0, 7.5, 3.2, 0.8, "Quality Filter\n& Dedup", COLORS["process"])
    draw_box(ax, 9.0, 7.5, 3.2, 0.8, "Candidate Site\nExtraction (BWA)", COLORS["process"])
    draw_box(ax, 15.0, 7.5, 3.2, 0.8, "Signal Normalisation\n(log1p + bin)", COLORS["process"])

    draw_arrow(ax, 3.0, 8.6, 3.0, 7.9)
    draw_arrow(ax, 9.0, 8.6, 9.0, 7.9)
    draw_arrow(ax, 15.0, 8.6, 15.0, 7.9)

    # ── Row 3: Feature Engineering ────────────────────────────────────────────
    draw_box(ax, 3.0, 6.0, 3.2, 0.8, "One-Hot Encoding\nGuide + Target (23×4)", COLORS["feature"])
    draw_box(ax, 9.0, 6.0, 3.2, 0.8, "Mismatch Pattern\n(23×15 channels)", COLORS["feature"])
    draw_box(ax, 15.0, 6.0, 3.2, 0.8, "Epigenetic Vector\n(8-dim)", COLORS["feature"])

    draw_arrow(ax, 3.0, 7.1, 3.0, 6.4)
    draw_arrow(ax, 9.0, 7.1, 9.0, 6.4)
    draw_arrow(ax, 15.0, 7.1, 15.0, 6.4)

    # ── Merge ─────────────────────────────────────────────────────────────────
    draw_box(ax, 9.0, 4.6, 4.0, 0.8, "Feature Fusion\n(23×23 seq tensor + 31-dim scalar)", COLORS["feature"])

    draw_arrow(ax, 3.0, 5.6, 6.5, 4.8)
    draw_arrow(ax, 9.0, 5.6, 9.0, 5.0)
    draw_arrow(ax, 15.0, 5.6, 11.5, 4.8)

    # ── Row 4: Model ──────────────────────────────────────────────────────────
    draw_box(ax, 3.5, 3.1, 3.0, 0.8, "Conv1D Block ×3\n(64→128→256 channels)", COLORS["model"])
    draw_box(ax, 7.5, 3.1, 2.8, 0.8, "Learnable\nPositional Encoding", COLORS["model"])
    draw_box(ax, 11.5, 3.1, 3.0, 0.8, "Multi-Head\nSelf-Attention (4 heads)", COLORS["model"])
    draw_box(ax, 15.5, 3.1, 2.5, 0.8, "Scalar\nEncoder (MLP)", COLORS["model"])

    draw_arrow(ax, 9.0, 4.2, 3.5, 3.5)
    draw_arrow(ax, 3.5, 2.7, 7.5, 3.1)
    draw_arrow(ax, 7.5, 2.7, 11.5, 3.1)
    draw_arrow(ax, 9.0, 4.2, 15.5, 3.5)

    # ── Row 5: Pooling + Fusion ───────────────────────────────────────────────
    draw_box(ax, 9.0, 1.8, 4.5, 0.8, "Global Avg+Max Pool → Concat with Scalar\n→ MLP Head → Sigmoid", COLORS["model"])
    draw_arrow(ax, 11.5, 2.7, 10.0, 2.2)
    draw_arrow(ax, 15.5, 2.7, 11.5, 2.0)

    # ── Row 6: Outputs ────────────────────────────────────────────────────────
    draw_box(ax, 4.5, 0.7, 2.8, 0.7, "Off-Target\nProbability", COLORS["output"])
    draw_box(ax, 9.0, 0.7, 2.8, 0.7, "AUROC / AUPRC\nEvaluation", COLORS["eval"])
    draw_box(ax, 13.5, 0.7, 2.8, 0.7, "SHAP Interpretability\n& Attention Maps", COLORS["eval"])

    draw_arrow(ax, 9.0, 1.4, 4.5, 1.05)
    draw_arrow(ax, 9.0, 1.4, 9.0, 1.05)
    draw_arrow(ax, 9.0, 1.4, 13.5, 1.05)

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_elements = [
        mpatches.Patch(facecolor=COLORS["input"],   label="Input Data"),
        mpatches.Patch(facecolor=COLORS["process"], label="Preprocessing"),
        mpatches.Patch(facecolor=COLORS["feature"], label="Feature Engineering"),
        mpatches.Patch(facecolor=COLORS["model"],   label="Model Components"),
        mpatches.Patch(facecolor=COLORS["output"],  label="Predictions"),
        mpatches.Patch(facecolor=COLORS["eval"],    label="Evaluation & Explainability"),
    ]
    ax.legend(handles=legend_elements, loc="lower left",
              fontsize=9, framealpha=0.8, ncol=3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Data flow diagram saved: {save_path}")


if __name__ == "__main__":
    create_dataflow_diagram()
