from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

FIGURE_SIZE = (16, 20)
PNG_OUTPUT = Path("/home/nahisaho/GitHub/aira/projects/ac30ac4e-cdba-463b-bde3-c9e7182c3e1b/workspace/figures/pipeline_architecture.png")
SVG_OUTPUT = Path("/home/nahisaho/GitHub/aira/projects/ac30ac4e-cdba-463b-bde3-c9e7182c3e1b/workspace/figures/pipeline_architecture.svg")
MANIFEST_OUTPUT = Path("/home/nahisaho/GitHub/aira/projects/ac30ac4e-cdba-463b-bde3-c9e7182c3e1b/workspace/results/pipeline_architecture_manifest.json")

MODULES = [
    (
        "Module 1: Signal-Level Basecalling (RNN)",
        [
            "• Bidirectional GRU layers",
            "• CTC decoding",
            "• Quality-aware consensus",
        ],
        "#DDEBF7",
    ),
    (
        "Module 2: Alignment & Feature Extraction",
        [
            "• Minimap2 alignment",
            "• Split-read signal extraction",
            "• Read-depth profiling",
            "• Soft-clip analysis",
        ],
        "#E7F4E4",
    ),
    (
        "Module 3: Integrated SV Detection",
        [
            "• Split-read caller",
            "• Read-depth caller",
            "• Local assembly caller (de novo)",
            "• Ensemble voting / evidence merging",
        ],
        "#FCE8D5",
    ),
    (
        "Module 4: Repeat Region Handler",
        [
            "• Telomere repeat detection (TTAGGG)",
            "• Centromere alpha-satellite analysis",
            "• Tandem repeat expansion detection",
            "• K-mer frequency filter",
        ],
        "#F9E2E7",
    ),
    (
        "Module 5: Complex SV Detector",
        [
            "• Chromothripsis pattern recognition",
            "• Extrachromosomal DNA (ecDNA) circular detection",
            "• Breakpoint graph construction",
            "• Multi-breakpoint clustering",
        ],
        "#E6E0F8",
    ),
    (
        "Module 6: Hybrid Integration",
        [
            "• Short-read evidence overlay",
            "• Genotype refinement",
            "• Breakpoint precision enhancement",
            "• Population frequency annotation",
        ],
        "#DFF2F1",
    ),
    (
        "Module 7: Quality & Output",
        [
            "• GIAB benchmark evaluation",
            "• Precision/Recall/F1 calculation",
            "• VCF output with confidence scores",
        ],
        "#F7F1D5",
    ),
]


def draw_box(ax, x, y, width, height, title, lines, facecolor, edgecolor="#4F6473"):
    shadow = FancyBboxPatch(
        (x + 0.008, y - 0.008),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=0,
        facecolor="#C7D2DA",
        alpha=0.35,
        zorder=1,
    )
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(shadow)
    ax.add_patch(box)

    ax.text(
        x + width / 2,
        y + height * 0.78,
        title,
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#24323D",
        zorder=3,
    )
    ax.text(
        x + 0.025,
        y + height * 0.47,
        "\n".join(lines),
        ha="left",
        va="center",
        fontsize=10.5,
        color="#2F3E46",
        zorder=3,
        linespacing=1.35,
    )


def draw_arrow(ax, start, end, color="#556B78", linestyle="solid", connectionstyle="arc3"):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            lw=2.0,
            color=color,
            linestyle=linestyle,
            shrinkA=6,
            shrinkB=6,
            mutation_scale=18,
            connectionstyle=connectionstyle,
        ),
        zorder=4,
    )


def build_figure():
    plt.rcParams["font.family"] = "DejaVu Sans"

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    fig.patch.set_facecolor("#F8FAFC")
    ax.set_facecolor("#F8FAFC")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.5,
        0.982,
        "DeepSV-LR: Long-Read Structural Variant Detection Pipeline",
        ha="center",
        va="top",
        fontsize=22,
        fontweight="bold",
        color="#1F2933",
    )
    ax.text(
        0.5,
        0.962,
        "Comprehensive architecture integrating signal, alignment, repeat-aware, complex SV, and hybrid evidence modules",
        ha="center",
        va="top",
        fontsize=11.5,
        color="#52606D",
    )
    ax.text(
        0.16,
        0.928,
        "Input Layer",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#334E68",
    )

    draw_box(
        ax,
        0.08,
        0.868,
        0.36,
        0.052,
        "Raw Signal",
        ["• ONT FAST5 / PacBio BAM"],
        "#E8EEF7",
    )
    draw_box(
        ax,
        0.56,
        0.868,
        0.30,
        0.052,
        "Short-read BAM",
        ["• Illumina"],
        "#F0E8F8",
    )

    x, width, height = 0.14, 0.72, 0.087
    start_y = 0.75
    gap = 0.105
    centers = []

    for idx, (title, lines, color) in enumerate(MODULES):
        y = start_y - idx * gap
        draw_box(ax, x, y, width, height, title, lines, color)
        centers.append((x + width / 2, y + height / 2, y, y + height))

    draw_arrow(ax, (0.26, 0.868), (centers[0][0], centers[0][3]))

    for idx in range(len(centers) - 1):
        draw_arrow(ax, (centers[idx][0], centers[idx][2]), (centers[idx + 1][0], centers[idx + 1][3]))

    draw_arrow(
        ax,
        (0.71, 0.868),
        (centers[5][0] + 0.12, centers[5][3]),
        color="#7B5EA7",
        linestyle="dashed",
        connectionstyle="angle3,angleA=-90,angleB=180",
    )

    ax.text(
        0.76,
        0.64,
        "Auxiliary\nshort-read evidence",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#6B4E8A",
        fontweight="bold",
    )

    ax.text(
        0.5,
        0.045,
        "Primary path: long-read signal to SV calling | Dashed path: hybrid evidence refinement",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#52606D",
    )

    return fig


def write_manifest():
    manifest = {
        "title": "DeepSV-LR pipeline architecture",
        "figure_size": list(FIGURE_SIZE),
        "png_output": str(PNG_OUTPUT),
        "svg_output": str(SVG_OUTPUT),
        "module_count": len(MODULES),
        "modules": [title for title, _, _ in MODULES],
        "inputs": ["Raw Signal (ONT FAST5 / PacBio BAM)", "Short-read BAM (Illumina)"],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main():
    PNG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    fig = build_figure()
    fig.savefig(PNG_OUTPUT, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(SVG_OUTPUT, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    write_manifest()
    print(f"Saved {PNG_OUTPUT}")
    print(f"Saved {SVG_OUTPUT}")
    print(f"Saved {MANIFEST_OUTPUT}")


if __name__ == "__main__":
    main()
