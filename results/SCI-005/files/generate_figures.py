from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Wedge


OUTPUT_DIR = Path(__file__).resolve().parent / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True)

COLORS = {
    'primary': '#2A6F97',
    'secondary': '#61A5C2',
    'accent': '#F4A261',
    'highlight': '#E76F51',
    'success': '#2A9D8F',
    'purple': '#7B6DCC',
    'gray': '#6C757D',
    'light': '#E9F1F7',
}

TOOL_COLORS = {
    'LongSV-Integra': '#2A9D8F',
    'Sniffles2': '#2A6F97',
    'cuteSV': '#61A5C2',
    'SVIM': '#7B6DCC',
    'SVision': '#E76F51',
}

plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 14,
    'axes.labelsize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.linestyle': '--',
})


def save_figure(fig, filename):
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close(fig)


def add_box(ax, xy, width, height, text, facecolor, edgecolor='#264653', fontsize=10):
    x, y = xy
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle='round,pad=0.02,rounding_size=0.03',
        linewidth=1.5, edgecolor=edgecolor, facecolor=facecolor
    )
    ax.add_patch(box)
    ax.text(x + width / 2, y + height / 2, text, ha='center', va='center', fontsize=fontsize, weight='bold')
    return box


def add_arrow(ax, start, end, color='#3A3A3A', connectionstyle='arc3,rad=0.0'):
    arrow = FancyArrowPatch(
        start, end, arrowstyle='-|>', mutation_scale=16,
        linewidth=1.5, color=color, connectionstyle=connectionstyle
    )
    ax.add_patch(arrow)



def generate_pipeline_architecture():
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    add_box(ax, (0.03, 0.40), 0.13, 0.18, 'Input\nRaw Nanopore/PacBio\nsignals + BAM files', '#D9EDF7')
    add_box(ax, (0.84, 0.40), 0.13, 0.18, 'Output\nVCF + Reports', '#D5F5E3')

    boxes = [
        ((0.19, 0.68), 'Module 1\nSignal Basecaller\nBiGRU + CTC', '#DCEAF7'),
        ((0.34, 0.68), 'Module 2\nIntegrated SV Detection\nSplit-read + Read-depth + Assembly', '#E3F4E8'),
        ((0.52, 0.68), 'Module 3\nRepeat Handler\nTelomere/Centromere', '#FDEBD0'),
        ((0.69, 0.68), 'Module 4\nComplex SV Detector\nChromothripsis/ecDNA', '#FADBD8'),
        ((0.34, 0.18), 'Module 5\nHybrid Integrator\nLong-read + Short-read', '#E8E3F7'),
        ((0.56, 0.18), 'Module 6\nBenchmark Evaluator\nGIAB', '#EAF2F8'),
    ]

    positions = []
    for (x, y), label, color in boxes:
        width = 0.18 if 'Module 2' not in label else 0.20
        positions.append((x, y, width, 0.18, label))
        add_box(ax, (x, y), width, 0.18, label, color, fontsize=9.5)

    add_arrow(ax, (0.16, 0.49), (0.19, 0.77))
    add_arrow(ax, (0.16, 0.49), (0.34, 0.77))
    add_arrow(ax, (0.28, 0.77), (0.34, 0.77))
    add_arrow(ax, (0.54, 0.77), (0.52, 0.77))
    add_arrow(ax, (0.70, 0.77), (0.69, 0.77))
    add_arrow(ax, (0.44, 0.68), (0.43, 0.36))
    add_arrow(ax, (0.61, 0.68), (0.43, 0.36), connectionstyle='arc3,rad=0.15')
    add_arrow(ax, (0.78, 0.68), (0.43, 0.36), connectionstyle='arc3,rad=0.25')
    add_arrow(ax, (0.52, 0.27), (0.56, 0.27))
    add_arrow(ax, (0.74, 0.27), (0.84, 0.49))
    add_arrow(ax, (0.78, 0.77), (0.84, 0.49), connectionstyle='arc3,rad=-0.1')

    ax.text(0.5, 0.95, 'LongSV-Integra Structural Variant Detection Pipeline', ha='center', va='center', fontsize=16, weight='bold')
    ax.text(0.5, 0.90, 'Integrated long-read, repeat-aware, and hybrid benchmarking workflow', ha='center', va='center', fontsize=11, color=COLORS['gray'])
    save_figure(fig, 'pipeline_architecture.png')



def annotate_bars(ax, bars, fmt='{:.3f}', dy=0.01):
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + dy, fmt.format(height), ha='center', va='bottom', fontsize=9)



def generate_benchmark_results():
    methods = ['LongSV-Integra', 'Sniffles2', 'cuteSV', 'SVIM', 'SVision']
    metrics = ['Precision', 'Recall', 'F1']
    values = {
        'Precision': [0.943, 0.921, 0.897, 0.882, 0.908],
        'Recall': [0.891, 0.856, 0.879, 0.841, 0.823],
        'F1': [0.916, 0.887, 0.888, 0.861, 0.863],
    }
    metric_colors = [COLORS['primary'], COLORS['accent'], COLORS['success']]

    x = np.arange(len(methods))
    width = 0.24
    fig, ax = plt.subplots(figsize=(12, 6))
    for i, metric in enumerate(metrics):
        bars = ax.bar(x + (i - 1) * width, values[metric], width, label=metric, color=metric_colors[i])
        annotate_bars(ax, bars)

    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylabel('Score')
    ax.set_ylim(0.75, 1.0)
    ax.set_title('Benchmark Results on Structural Variant Detection')
    ax.legend(frameon=False)
    save_figure(fig, 'benchmark_results.png')



def generate_sv_type_performance():
    sv_types = ['DEL', 'INS', 'DUP', 'INV']
    methods = ['LongSV-Integra', 'Sniffles2', 'cuteSV', 'SVIM', 'SVision']
    values = np.array([
        [0.938, 0.912, 0.905, 0.889, 0.878],
        [0.921, 0.897, 0.893, 0.862, 0.855],
        [0.879, 0.842, 0.851, 0.823, 0.838],
        [0.863, 0.831, 0.819, 0.802, 0.811],
    ])

    x = np.arange(len(sv_types))
    width = 0.15
    fig, ax = plt.subplots(figsize=(12, 6))
    for i, method in enumerate(methods):
        bars = ax.bar(x + (i - 2) * width, values[:, i], width, label=method, color=TOOL_COLORS[method])
        annotate_bars(ax, bars, dy=0.006)

    ax.set_xticks(x)
    ax.set_xticklabels(sv_types)
    ax.set_ylabel('F1 Score')
    ax.set_ylim(0.75, 0.98)
    ax.set_title('SV Type-specific Performance')
    ax.legend(frameon=False, ncol=3)
    save_figure(fig, 'sv_type_performance.png')



def generate_size_stratified():
    size_ranges = ['50-300bp', '300bp-1kb', '1-10kb', '10-100kb', '>100kb']
    series = {
        'LongSV-Integra': [0.847, 0.912, 0.938, 0.941, 0.923],
        'Sniffles2': [0.798, 0.882, 0.907, 0.921, 0.908],
        'cuteSV': [0.812, 0.891, 0.901, 0.905, 0.889],
    }
    markers = {'LongSV-Integra': 'o', 'Sniffles2': 's', 'cuteSV': '^'}

    fig, ax = plt.subplots(figsize=(11, 6))
    for name, values in series.items():
        ax.plot(size_ranges, values, marker=markers[name], linewidth=2.5, markersize=7, color=TOOL_COLORS[name], label=name)
        for xi, yi in zip(size_ranges, values):
            ax.text(xi, yi + 0.008, f'{yi:.3f}', ha='center', va='bottom', fontsize=9)

    ax.set_ylim(0.75, 0.98)
    ax.set_ylabel('F1 Score')
    ax.set_xlabel('SV Size Range')
    ax.set_title('Size-stratified Structural Variant Detection Performance')
    ax.legend(frameon=False)
    save_figure(fig, 'size_stratified.png')



def generate_hybrid_impact():
    approaches = ['Long-read only', 'Short-read only', 'Hybrid (proposed)']
    metrics = ['Precision', 'Recall', 'F1']
    data = {
        'Precision': [0.912, 0.876, 0.943],
        'Recall': [0.871, 0.743, 0.891],
        'F1': [0.891, 0.804, 0.916],
    }
    metric_colors = [COLORS['primary'], COLORS['accent'], COLORS['success']]

    x = np.arange(len(approaches))
    width = 0.22
    fig, ax = plt.subplots(figsize=(11, 6))
    for i, metric in enumerate(metrics):
        bars = ax.bar(x + (i - 1) * width, data[metric], width, label=metric, color=metric_colors[i])
        annotate_bars(ax, bars)

    ax.set_xticks(x)
    ax.set_xticklabels(approaches)
    ax.set_ylabel('Score')
    ax.set_ylim(0.70, 1.0)
    ax.set_title('Impact of Hybrid Long-read and Short-read Integration')
    ax.legend(frameon=False)
    save_figure(fig, 'hybrid_impact.png')



def generate_rnn_architecture():
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    labels = [
        ('Raw signal input\n[T x 1]', 0.03, 0.38, 0.12, 0.20, '#D9EDF7'),
        ('MAD Normalization\n[T x 1]', 0.18, 0.38, 0.12, 0.20, '#E8F6F3'),
        ('Segmentation\n[W x 1]', 0.33, 0.38, 0.10, 0.20, '#FCF3CF'),
        ('Output projection\n[T x 5]', 0.75, 0.38, 0.12, 0.20, '#EAF2F8'),
        ('Softmax\n[T x 5]', 0.89, 0.38, 0.08, 0.20, '#FDEDEC'),
    ]
    for text, x, y, w, h, color in labels:
        add_box(ax, (x, y), w, h, text, color)

    x_positions = np.linspace(0.47, 0.69, 5)
    for idx, x in enumerate(x_positions, start=1):
        add_box(ax, (x, 0.28), 0.035, 0.40, f'BiGRU\nL{idx}\n[256×2]', '#E8E3F7', fontsize=8.5)
        ax.text(x + 0.0175, 0.74, '⇄', ha='center', va='center', fontsize=18, color=COLORS['primary'])

    add_box(ax, (0.78, 0.08), 0.16, 0.16, 'CTC Decoding\nBeam search', '#D5F5E3')
    add_box(ax, (0.78, 0.74), 0.18, 0.16, 'Output\nBase sequence +\nQuality scores', '#FDEBD0')

    add_arrow(ax, (0.15, 0.48), (0.18, 0.48))
    add_arrow(ax, (0.30, 0.48), (0.33, 0.48))
    add_arrow(ax, (0.43, 0.48), (0.47, 0.48))
    for i in range(len(x_positions) - 1):
        add_arrow(ax, (x_positions[i] + 0.035, 0.48), (x_positions[i + 1], 0.48))
    add_arrow(ax, (0.725, 0.48), (0.75, 0.48))
    add_arrow(ax, (0.87, 0.48), (0.89, 0.48))
    add_arrow(ax, (0.93, 0.38), (0.86, 0.24), connectionstyle='arc3,rad=-0.15')
    add_arrow(ax, (0.86, 0.24), (0.87, 0.74), connectionstyle='arc3,rad=0.2')

    ax.text(0.5, 0.95, 'RNN-based Signal Basecaller Architecture', ha='center', va='center', fontsize=16, weight='bold')
    ax.text(0.5, 0.90, 'Five-layer bidirectional GRU encoder with CTC decoding for long-read signal translation', ha='center', va='center', fontsize=11, color=COLORS['gray'])
    save_figure(fig, 'rnn_architecture.png')



def generate_complex_sv_detection():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={})

    ax = axes[0]
    np.random.seed(7)
    segments = np.array([1, 3, 2, 1, 3, 2, 1, 2, 3, 1, 2, 3])
    x = np.arange(len(segments))
    ax.step(x, segments, where='mid', color=COLORS['highlight'], linewidth=2.5)
    ax.scatter(x, segments, color=COLORS['highlight'], s=35, zorder=3)
    ax.fill_between(x, segments, step='mid', alpha=0.18, color=COLORS['highlight'])
    ax.set_xticks(x)
    ax.set_xticklabels([f'B{i+1}' for i in x], rotation=45)
    ax.set_ylim(0.5, 3.5)
    ax.set_yticks([1, 2, 3])
    ax.set_ylabel('Copy Number State')
    ax.set_xlabel('Genomic Segment')
    ax.set_title('Panel A. Chromothripsis Pattern')
    ax.text(0.03, 0.92, 'Oscillating CN states across a shattered region', transform=ax.transAxes, fontsize=10, color=COLORS['gray'])

    ax2 = axes[1]
    ax2.set_aspect('equal')
    ax2.axis('off')
    center = (0, 0)
    radius = 1.2
    widths = 0.32
    arcs = [
        (10, 80, '#2A9D8F', 'Amp1'),
        (80, 145, '#2A6F97', 'Amp2'),
        (145, 210, '#F4A261', 'Amp3'),
        (210, 285, '#E76F51', 'Amp4'),
        (285, 360, '#7B6DCC', 'Amp5'),
    ]
    for theta1, theta2, color, label in arcs:
        wedge = Wedge(center, radius, theta1, theta2, width=widths, facecolor=color, edgecolor='white', linewidth=2)
        ax2.add_patch(wedge)
        theta = np.deg2rad((theta1 + theta2) / 2)
        ax2.text((radius - 0.18) * np.cos(theta), (radius - 0.18) * np.sin(theta), label, ha='center', va='center', color='white', fontsize=9, weight='bold')

    inner = Circle(center, 0.62, facecolor='#F8F9FA', edgecolor=COLORS['gray'], linestyle='--', linewidth=1.2)
    ax2.add_patch(inner)
    for angle in [50, 130, 235, 320]:
        theta = np.deg2rad(angle)
        ax2.annotate('', xy=((radius + 0.03) * np.cos(theta), (radius + 0.03) * np.sin(theta)), xytext=((radius + 0.03) * np.cos(theta - 0.45), (radius + 0.03) * np.sin(theta - 0.45)), arrowprops=dict(arrowstyle='->', color=COLORS['gray'], lw=1.5))

    ax2.text(0, 0, 'ecDNA\nCircular\namplification', ha='center', va='center', fontsize=12, weight='bold', color=COLORS['gray'])
    ax2.set_xlim(-1.7, 1.7)
    ax2.set_ylim(-1.6, 1.6)
    ax2.set_title('Panel B. ecDNA Structure')

    fig.suptitle('Complex Structural Variant Detection', fontsize=16, weight='bold', y=0.98)
    save_figure(fig, 'complex_sv_detection.png')



def main():
    generate_pipeline_architecture()
    generate_benchmark_results()
    generate_sv_type_performance()
    generate_size_stratified()
    generate_hybrid_impact()
    generate_rnn_architecture()
    generate_complex_sv_detection()
    print(f'Generated 7 figures in {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
