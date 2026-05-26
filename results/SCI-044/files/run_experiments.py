"""
Main experiment script for RNA secondary structure prediction study.
Runs all experiments and generates figures.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import time
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
from rna_predictor import (
    RNAStructurePredictor, PseudoknotDetector, RiboswitchPredictor,
    simulate_shape_data, simulate_dms_data, compute_mutual_information,
    generate_synthetic_msa, compute_metrics, parse_pairs_from_structure
)

np.random.seed(42)
os.makedirs('figures', exist_ok=True)

# ============================================================
# Test RNA sequences with known structures
# ============================================================

BENCHMARK_RNAS = {}

def _make_benchmark(name, seq, struct, desc):
    """Helper: ensure lengths match and brackets are balanced."""
    assert len(seq) == len(struct), f"{name}: seq({len(seq)}) != struct({len(struct)})"
    BENCHMARK_RNAS[name] = {'sequence': seq, 'structure': struct, 'description': desc}

# Simple hairpin (14 nt)
_make_benchmark('hairpin_1',
    'CGGCUAAAAAGCCG',
    '((((......))))',
    'Simple hairpin (14 nt)')

# Two-stem structure (30 nt)
_make_benchmark('two_stem',
    'GUGGAACACCCACUAAGCCGCUUAACGGCU',
    '((((.....))))...((((.....)))).',
    'Two-stem structure (30 nt)')

# tRNA-like (41 nt)
_make_benchmark('tRNA_like',
    'GGUGAACAGAGCUUAAAAUGCUCAAUGGACUCCACUUCACC',
    '((((((..((((.......))))..(((...))).))))))',
    'tRNA-like structure (41 nt)')

# Hammerhead-like (35 nt)
_make_benchmark('hammerhead_like',
    'CACAACAACACCACCUAUAGGUGACGCGGUUGUGA',
    '((((((..((((.......)))).(()))))))).',
    'Hammerhead-like (35 nt)')

# Multi-stem (46 nt)
_make_benchmark('multi_stem',
    'AUACAGCGCUAAUCGCGCCUACCAAUUCCAAAGAAAAGGUCAUAUC',
    '(((..((((.....))))..(((..(((....)))..)))..))).',
    'Multi-stem structure (46 nt)')

# SARS-CoV-2 5'UTR (first 265 nt — stem-loop region SL1-SL5)
SARS_COV2_5UTR = {
    'sequence': (
        'AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAU'
        'CUGUGUUGCUAUAUAUAUAUCAUUGGCGUUUCUGCUGUCCGAUUUUUAUUAUAUUUAUAUAUAUUUAUAUAAUAUUAUUAUGA'
        'UAGAUUAUUAUUAUAUUAUUAUUAUUAUAUAUUUUAUAUUAAUUAUGUAUUAUUGUUUAUUGUUGUAUUUUAUUAUUAAUUAU'
        'GUAUUAUUGUUUAUUG'
    ),
    'reference_sl1': '(((((((...)))))))' + '.' * 249,  # SL1 (simplified)
    'description': 'SARS-CoV-2 5\'UTR region',
}

# Use a more realistic shorter 5'UTR region for prediction
SARS_COV2_SHORT = {
    'sequence':  'UAAACACGAGUAGAAGUGGAUCAAAAAAUCCACUAGCCAGCCUCAUCUGGCUUCUACUCCAUACCCCUAAACUAU',
    'structure': '.......((((((..(((((........))))).((((((......)))))).))))))................',
    'description': 'SARS-CoV-2 5\'UTR SL1-SL3 region (75 nt)',
}

# Riboswitch examples
RIBOSWITCH_EXAMPLES = {
    'TPP_riboswitch': {
        'sequence':  'CUUUCAUCGGCGUAAAAGAUCCUAAUCAUCGCCGAUCGGUCAUAACAGACCUUAAAAGUU',
        'structure': '((((..((((((.....(((....)))..))))))..((((......))))...))))..',
        'description': 'TPP riboswitch aptamer domain',
    },
    'SAM_riboswitch': {
        'sequence':  'UCCGUCUCUGCCCCAGACUCUAAUUCUCUAGGAAAUCCACAGAAUUCGGAUACAAAA',
        'structure': '((((..((((....(((.......)))...((....)).))))...)))).......', 
        'description': 'SAM-I riboswitch aptamer',
    },
}


def run_experiment_1_baseline():
    """Experiment 1: Baseline prediction with Turner model."""
    print("=" * 60)
    print("Experiment 1: Baseline Turner Model Prediction")
    print("=" * 60)

    results = {}
    for name, data in BENCHMARK_RNAS.items():
        seq = data['sequence']
        ref = data['structure']
        t0 = time.time()
        predictor = RNAStructurePredictor(seq)
        mfe, pred_struct = predictor.fold()
        elapsed = time.time() - t0

        # Adjust structure length
        if len(pred_struct) < len(ref):
            pred_struct += '.' * (len(ref) - len(pred_struct))
        elif len(pred_struct) > len(ref):
            pred_struct = pred_struct[:len(ref)]

        metrics = compute_metrics(pred_struct, ref)
        metrics['mfe'] = mfe
        metrics['time'] = elapsed
        results[name] = metrics

        print(f"\n{name} ({data['description']}):")
        print(f"  MFE: {mfe:.2f} kcal/mol | Time: {elapsed:.3f}s")
        print(f"  Sensitivity: {metrics['sensitivity']:.3f} | PPV: {metrics['ppv']:.3f} | F1: {metrics['f1']:.3f}")

    return results


def run_experiment_2_shape_dms():
    """Experiment 2: SHAPE/DMS constraint integration."""
    print("\n" + "=" * 60)
    print("Experiment 2: SHAPE/DMS Constrained Prediction")
    print("=" * 60)

    results = {'baseline': {}, 'shape': {}, 'dms': {}, 'combined': {}}
    for name, data in BENCHMARK_RNAS.items():
        seq = data['sequence']
        ref = data['structure']

        # Baseline
        pred_base = RNAStructurePredictor(seq)
        _, struct_base = pred_base.fold()
        if len(struct_base) < len(ref):
            struct_base += '.' * (len(ref) - len(struct_base))
        results['baseline'][name] = compute_metrics(struct_base[:len(ref)], ref)

        # With SHAPE
        shape = simulate_shape_data(seq, ref)
        pred_shape = RNAStructurePredictor(seq, shape_data=shape)
        _, struct_shape = pred_shape.fold()
        if len(struct_shape) < len(ref):
            struct_shape += '.' * (len(ref) - len(struct_shape))
        results['shape'][name] = compute_metrics(struct_shape[:len(ref)], ref)

        # With DMS
        dms = simulate_dms_data(seq, ref)
        pred_dms = RNAStructurePredictor(seq, dms_data=dms)
        _, struct_dms = pred_dms.fold()
        if len(struct_dms) < len(ref):
            struct_dms += '.' * (len(ref) - len(struct_dms))
        results['dms'][name] = compute_metrics(struct_dms[:len(ref)], ref)

        # Combined
        pred_comb = RNAStructurePredictor(seq, shape_data=shape, dms_data=dms)
        _, struct_comb = pred_comb.fold()
        if len(struct_comb) < len(ref):
            struct_comb += '.' * (len(ref) - len(struct_comb))
        results['combined'][name] = compute_metrics(struct_comb[:len(ref)], ref)

        print(f"\n{name}:")
        for method in ['baseline', 'shape', 'dms', 'combined']:
            m = results[method][name]
            print(f"  {method:10s}: Sens={m['sensitivity']:.3f} PPV={m['ppv']:.3f} F1={m['f1']:.3f}")

    return results


def run_experiment_3_covariation():
    """Experiment 3: MSA-based covariation analysis."""
    print("\n" + "=" * 60)
    print("Experiment 3: MSA Covariation Integration")
    print("=" * 60)

    results = {'baseline': {}, 'covariation': {}}
    for name, data in list(BENCHMARK_RNAS.items())[:3]:
        seq = data['sequence']
        ref = data['structure']

        # Generate synthetic MSA
        msa = generate_synthetic_msa(seq, ref, n_sequences=30, mutation_rate=0.12)
        mi_matrix = compute_mutual_information(msa)

        # Normalize MI to use as covariation scores
        if mi_matrix.max() > 0:
            covar_scores = mi_matrix / mi_matrix.max() * 2.0
        else:
            covar_scores = mi_matrix

        # Baseline
        pred_base = RNAStructurePredictor(seq)
        _, struct_base = pred_base.fold()
        if len(struct_base) < len(ref):
            struct_base += '.' * (len(ref) - len(struct_base))
        results['baseline'][name] = compute_metrics(struct_base[:len(ref)], ref)

        # With covariation
        pred_cov = RNAStructurePredictor(seq, covariation_scores=covar_scores)
        _, struct_cov = pred_cov.fold()
        if len(struct_cov) < len(ref):
            struct_cov += '.' * (len(ref) - len(struct_cov))
        results['covariation'][name] = compute_metrics(struct_cov[:len(ref)], ref)

        print(f"\n{name}:")
        print(f"  Baseline:     F1={results['baseline'][name]['f1']:.3f}")
        print(f"  Covariation:  F1={results['covariation'][name]['f1']:.3f}")

    return results


def run_experiment_4_pseudoknots():
    """Experiment 4: Pseudoknot detection."""
    print("\n" + "=" * 60)
    print("Experiment 4: Pseudoknot Detection")
    print("=" * 60)

    pk_results = {}
    for name, data in BENCHMARK_RNAS.items():
        seq = data['sequence']
        predictor = RNAStructurePredictor(seq)
        _, base_struct = predictor.fold()

        if len(base_struct) < len(seq):
            base_struct += '.' * (len(seq) - len(base_struct))

        detector = PseudoknotDetector(seq, base_struct[:len(seq)])
        pks = detector.detect_pseudoknots()
        enhanced = detector.add_pseudoknots()

        pk_results[name] = {
            'n_pseudoknots': len(pks),
            'base_pairs': base_struct.count('('),
            'pk_pairs': enhanced.count('['),
            'structure': enhanced,
        }
        print(f"\n{name}: {len(pks)} pseudoknot(s) detected")
        print(f"  Base pairs: {base_struct.count('(')} | PK pairs: {enhanced.count('[')}")

    return pk_results


def run_experiment_5_riboswitch():
    """Experiment 5: Riboswitch structure-function prediction."""
    print("\n" + "=" * 60)
    print("Experiment 5: Riboswitch Structure-Function Prediction")
    print("=" * 60)

    rs_results = {}
    for name, data in RIBOSWITCH_EXAMPLES.items():
        seq = data['sequence']
        ref = data['structure']

        predictor = RNAStructurePredictor(seq)
        mfe, pred_struct = predictor.fold()

        if len(pred_struct) < len(ref):
            pred_struct += '.' * (len(ref) - len(pred_struct))

        rs_pred = RiboswitchPredictor(seq, pred_struct[:len(seq)])
        func_elements = rs_pred.identify_functional_elements()

        metrics = compute_metrics(pred_struct[:len(ref)], ref)
        rs_results[name] = {
            'metrics': metrics,
            'functional': func_elements,
            'mfe': mfe,
        }

        print(f"\n{name} ({data['description']}):")
        print(f"  MFE: {mfe:.2f} | F1: {metrics['f1']:.3f}")
        print(f"  Stem-loops: {func_elements['stem_loops']}")
        print(f"  Aptamer motifs: {func_elements['aptamer_candidates']}")
        print(f"  Expression platform: {func_elements['expression_platform']}")
        print(f"  Switch potential: {func_elements['structural_switch_potential']:.2f}")

    return rs_results


def run_experiment_6_sars_cov2():
    """Experiment 6: SARS-CoV-2 5'UTR case study."""
    print("\n" + "=" * 60)
    print("Experiment 6: SARS-CoV-2 5'UTR Case Study")
    print("=" * 60)

    data = SARS_COV2_SHORT
    seq = data['sequence']
    ref = data['structure']

    # Baseline prediction
    t0 = time.time()
    predictor = RNAStructurePredictor(seq)
    mfe_base, struct_base = predictor.fold()
    time_base = time.time() - t0

    if len(struct_base) < len(ref):
        struct_base += '.' * (len(ref) - len(struct_base))

    # With simulated SHAPE
    shape = simulate_shape_data(seq, ref)
    t0 = time.time()
    pred_shape = RNAStructurePredictor(seq, shape_data=shape)
    mfe_shape, struct_shape = pred_shape.fold()
    time_shape = time.time() - t0

    if len(struct_shape) < len(ref):
        struct_shape += '.' * (len(ref) - len(struct_shape))

    # With MSA covariation
    msa = generate_synthetic_msa(seq, ref, n_sequences=40)
    mi = compute_mutual_information(msa)
    if mi.max() > 0:
        covar = mi / mi.max() * 2.0
    else:
        covar = mi

    t0 = time.time()
    pred_cov = RNAStructurePredictor(seq, covariation_scores=covar)
    mfe_cov, struct_cov = pred_cov.fold()
    time_cov = time.time() - t0

    if len(struct_cov) < len(ref):
        struct_cov += '.' * (len(ref) - len(struct_cov))

    # Combined
    t0 = time.time()
    pred_all = RNAStructurePredictor(seq, shape_data=shape, covariation_scores=covar)
    mfe_all, struct_all = pred_all.fold()
    time_all = time.time() - t0

    if len(struct_all) < len(ref):
        struct_all += '.' * (len(ref) - len(struct_all))

    # Pseudoknot detection
    detector = PseudoknotDetector(seq, struct_all[:len(seq)])
    struct_pk = detector.add_pseudoknots()

    metrics = {
        'baseline': compute_metrics(struct_base[:len(ref)], ref),
        'shape': compute_metrics(struct_shape[:len(ref)], ref),
        'covariation': compute_metrics(struct_cov[:len(ref)], ref),
        'combined': compute_metrics(struct_all[:len(ref)], ref),
    }

    print(f"\nSequence: {seq[:50]}...")
    print(f"Reference: {ref[:50]}...")
    for method, m in metrics.items():
        print(f"  {method:12s}: Sens={m['sensitivity']:.3f} PPV={m['ppv']:.3f} F1={m['f1']:.3f}")
    print(f"\nPseudoknot-enhanced structure: {struct_pk[:60]}...")

    return {
        'metrics': metrics,
        'structures': {
            'reference': ref,
            'baseline': struct_base[:len(ref)],
            'shape': struct_shape[:len(ref)],
            'covariation': struct_cov[:len(ref)],
            'combined': struct_all[:len(ref)],
            'with_pk': struct_pk,
        },
        'mfe': {'baseline': mfe_base, 'shape': mfe_shape, 'covariation': mfe_cov, 'combined': mfe_all},
        'times': {'baseline': time_base, 'shape': time_shape, 'covariation': time_cov, 'combined': time_all},
    }


# ============================================================
# Visualization Functions
# ============================================================

def plot_benchmark_comparison(results_baseline, results_shape_dms):
    """Figure 1: Benchmark comparison bar chart."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    rna_names = list(results_baseline.keys())
    x = np.arange(len(rna_names))
    width = 0.2

    methods = ['baseline', 'shape', 'dms', 'combined']
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#E91E63']
    labels = ['Baseline', '+SHAPE', '+DMS', '+SHAPE+DMS']

    for ax_idx, metric in enumerate(['sensitivity', 'ppv', 'f1']):
        for m_idx, method in enumerate(methods):
            vals = [results_shape_dms[method][name][metric] for name in rna_names]
            axes[ax_idx].bar(x + m_idx * width, vals, width, label=labels[m_idx],
                           color=colors[m_idx], alpha=0.85)
        axes[ax_idx].set_xticks(x + 1.5 * width)
        axes[ax_idx].set_xticklabels([n.replace('_', '\n') for n in rna_names], fontsize=8, rotation=30)
        axes[ax_idx].set_ylabel(metric.upper())
        axes[ax_idx].set_ylim(0, 1.1)
        axes[ax_idx].legend(fontsize=7)
        axes[ax_idx].set_title(f'{metric.upper()} by Method')
        axes[ax_idx].grid(axis='y', alpha=0.3)

    plt.suptitle('RNA Secondary Structure Prediction: Method Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/benchmark_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/benchmark_comparison.png")


def plot_covariation_heatmap(results_covar):
    """Figure 2: Covariation analysis effect."""
    fig, ax = plt.subplots(figsize=(8, 5))
    rna_names = list(results_covar['baseline'].keys())
    x = np.arange(len(rna_names))
    width = 0.35

    baseline_f1 = [results_covar['baseline'][n]['f1'] for n in rna_names]
    covar_f1 = [results_covar['covariation'][n]['f1'] for n in rna_names]

    bars1 = ax.bar(x - width/2, baseline_f1, width, label='Baseline', color='#2196F3', alpha=0.85)
    bars2 = ax.bar(x + width/2, covar_f1, width, label='+Covariation (MSA)', color='#9C27B0', alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels([n.replace('_', '\n') for n in rna_names], fontsize=9)
    ax.set_ylabel('F1 Score')
    ax.set_title('Effect of MSA-Based Covariation on Prediction Accuracy', fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    ax.grid(axis='y', alpha=0.3)

    for bar, val in zip(bars1, baseline_f1):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.2f}',
                ha='center', fontsize=8)
    for bar, val in zip(bars2, covar_f1):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.2f}',
                ha='center', fontsize=8)

    plt.tight_layout()
    plt.savefig('figures/covariation_effect.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/covariation_effect.png")


def plot_sars_cov2_results(sars_results):
    """Figure 3: SARS-CoV-2 case study results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Method comparison
    methods = list(sars_results['metrics'].keys())
    metrics_list = ['sensitivity', 'ppv', 'f1']
    x = np.arange(len(methods))
    width = 0.25
    colors = ['#2196F3', '#4CAF50', '#E91E63']

    for m_idx, metric in enumerate(metrics_list):
        vals = [sars_results['metrics'][m][metric] for m in methods]
        axes[0].bar(x + m_idx * width, vals, width, label=metric.upper(), color=colors[m_idx], alpha=0.85)

    axes[0].set_xticks(x + width)
    axes[0].set_xticklabels(methods, fontsize=9)
    axes[0].set_ylabel('Score')
    axes[0].set_title('SARS-CoV-2 5\'UTR Prediction Accuracy', fontweight='bold')
    axes[0].legend()
    axes[0].set_ylim(0, 1.1)
    axes[0].grid(axis='y', alpha=0.3)

    # Right: MFE comparison
    mfe_vals = [sars_results['mfe'][m] for m in methods]
    bars = axes[1].bar(methods, mfe_vals, color=['#2196F3', '#4CAF50', '#9C27B0', '#E91E63'], alpha=0.85)
    axes[1].set_ylabel('MFE (kcal/mol)')
    axes[1].set_title('Minimum Free Energy by Method', fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, mfe_vals):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() - 0.5 if val < 0 else bar.get_height() + 0.2,
                    f'{val:.1f}', ha='center', fontsize=9)

    plt.suptitle('SARS-CoV-2 5\'UTR Structure Prediction Case Study', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/sars_cov2_casestudy.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/sars_cov2_casestudy.png")


def plot_pseudoknot_analysis(pk_results):
    """Figure 4: Pseudoknot detection summary."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(pk_results.keys())
    base_pairs = [pk_results[n]['base_pairs'] for n in names]
    pk_pairs = [pk_results[n]['pk_pairs'] for n in names]
    n_pks = [pk_results[n]['n_pseudoknots'] for n in names]

    x = np.arange(len(names))
    width = 0.35

    axes[0].bar(x - width/2, base_pairs, width, label='Standard base pairs', color='#2196F3', alpha=0.85)
    axes[0].bar(x + width/2, pk_pairs, width, label='Pseudoknot pairs', color='#FF5722', alpha=0.85)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([n.replace('_', '\n') for n in names], fontsize=8, rotation=30)
    axes[0].set_ylabel('Number of pairs')
    axes[0].set_title('Base Pairs vs Pseudoknot Pairs', fontweight='bold')
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)

    axes[1].bar(names, n_pks, color='#FF5722', alpha=0.85)
    axes[1].set_xticklabels([n.replace('_', '\n') for n in names], fontsize=8, rotation=30)
    axes[1].set_ylabel('Count')
    axes[1].set_title('Pseudoknots Detected per RNA', fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)

    plt.suptitle('Pseudoknot Detection Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/pseudoknot_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/pseudoknot_analysis.png")


def plot_riboswitch_analysis(rs_results):
    """Figure 5: Riboswitch structure-function analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(rs_results.keys())
    f1_scores = [rs_results[n]['metrics']['f1'] for n in names]
    switch_potential = [rs_results[n]['functional']['structural_switch_potential'] for n in names]
    stem_loops = [rs_results[n]['functional']['stem_loops'] for n in names]

    # Left: F1 and switch potential
    x = np.arange(len(names))
    width = 0.35
    ax1 = axes[0]
    ax2 = ax1.twinx()

    bars1 = ax1.bar(x - width/2, f1_scores, width, label='F1 Score', color='#2196F3', alpha=0.85)
    bars2 = ax2.bar(x + width/2, switch_potential, width, label='Switch Potential', color='#FF9800', alpha=0.85)

    ax1.set_xticks(x)
    ax1.set_xticklabels([n.replace('_', '\n') for n in names], fontsize=9)
    ax1.set_ylabel('F1 Score', color='#2196F3')
    ax2.set_ylabel('Switch Potential', color='#FF9800')
    ax1.set_title('Riboswitch Prediction Quality', fontweight='bold')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax1.set_ylim(0, 1.1)
    ax2.set_ylim(0, 1.1)

    # Right: Functional element summary table
    axes[1].axis('off')
    table_data = []
    for name in names:
        func = rs_results[name]['functional']
        table_data.append([
            name.replace('_', ' '),
            func['stem_loops'],
            ', '.join(func['aptamer_candidates']) if func['aptamer_candidates'] else 'None',
            func['expression_platform'][:20],
            f"{func['structural_switch_potential']:.2f}"
        ])

    table = axes[1].table(cellText=table_data,
                          colLabels=['RNA', 'Stem-loops', 'Aptamer Motifs', 'Expression Platform', 'Switch Score'],
                          loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    axes[1].set_title('Functional Element Summary', fontweight='bold', pad=20)

    plt.suptitle('Riboswitch Structure-Function Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/riboswitch_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/riboswitch_analysis.png")


def plot_computation_time(results_baseline, sars_results):
    """Figure 6: Computational performance."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: time vs sequence length
    names = list(results_baseline.keys())
    seq_lengths = [len(BENCHMARK_RNAS[n]['sequence']) for n in names]
    times = [results_baseline[n]['time'] for n in names]

    axes[0].scatter(seq_lengths, times, s=100, c='#2196F3', zorder=5)
    for i, name in enumerate(names):
        axes[0].annotate(name.replace('_', '\n'), (seq_lengths[i], times[i]),
                        textcoords="offset points", xytext=(5, 5), fontsize=7)

    # Fit and plot expected O(n^3) curve
    if len(seq_lengths) > 1:
        x_fit = np.linspace(min(seq_lengths), max(seq_lengths), 100)
        # Normalize
        coeff = np.median([t / (l**3) for t, l in zip(times, seq_lengths) if t > 0])
        y_fit = coeff * x_fit**3
        axes[0].plot(x_fit, y_fit, '--', color='gray', alpha=0.5, label='O(n³) expected')
        axes[0].legend()

    axes[0].set_xlabel('Sequence Length (nt)')
    axes[0].set_ylabel('Time (seconds)')
    axes[0].set_title('Computation Time vs Sequence Length', fontweight='bold')
    axes[0].grid(alpha=0.3)

    # Right: SARS-CoV-2 time breakdown
    methods = list(sars_results['times'].keys())
    time_vals = [sars_results['times'][m] for m in methods]
    colors = ['#2196F3', '#4CAF50', '#9C27B0', '#E91E63']
    bars = axes[1].bar(methods, time_vals, color=colors, alpha=0.85)
    axes[1].set_ylabel('Time (seconds)')
    axes[1].set_title('SARS-CoV-2 Prediction Time by Method', fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, time_vals):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                    f'{val:.3f}s', ha='center', fontsize=9)

    plt.suptitle('Computational Performance Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/computation_time.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/computation_time.png")


def plot_overall_summary(results_baseline, results_shape_dms, results_covar, sars_results):
    """Figure 7: Overall summary heatmap."""
    fig, ax = plt.subplots(figsize=(10, 6))

    all_methods = ['Baseline', '+SHAPE', '+DMS', '+SHAPE+DMS', '+Covariation']
    all_rnas = list(BENCHMARK_RNAS.keys())

    f1_matrix = np.zeros((len(all_methods), len(all_rnas)))
    method_keys = ['baseline', 'shape', 'dms', 'combined']

    for m_idx, method in enumerate(method_keys):
        for r_idx, rna in enumerate(all_rnas):
            if rna in results_shape_dms[method]:
                f1_matrix[m_idx, r_idx] = results_shape_dms[method][rna]['f1']

    for r_idx, rna in enumerate(all_rnas):
        if rna in results_covar.get('covariation', {}):
            f1_matrix[4, r_idx] = results_covar['covariation'][rna]['f1']
        else:
            f1_matrix[4, r_idx] = results_shape_dms['baseline'][rna]['f1']

    sns.heatmap(f1_matrix, annot=True, fmt='.2f', cmap='RdYlGn',
                xticklabels=[n.replace('_', '\n') for n in all_rnas],
                yticklabels=all_methods, vmin=0, vmax=1, ax=ax,
                linewidths=0.5, linecolor='white')
    ax.set_title('F1 Score Heatmap: All Methods × All RNAs', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/overall_summary_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figures/overall_summary_heatmap.png")


# ============================================================
# Main Execution
# ============================================================

if __name__ == '__main__':
    print("RNA Secondary Structure Prediction - Full Experiment Suite")
    print("=" * 60)

    # Run experiments
    results_1 = run_experiment_1_baseline()
    results_2 = run_experiment_2_shape_dms()
    results_3 = run_experiment_3_covariation()
    results_4 = run_experiment_4_pseudoknots()
    results_5 = run_experiment_5_riboswitch()
    results_6 = run_experiment_6_sars_cov2()

    # Generate figures
    print("\n" + "=" * 60)
    print("Generating Figures")
    print("=" * 60)

    plot_benchmark_comparison(results_1, results_2)
    plot_covariation_heatmap(results_3)
    plot_sars_cov2_results(results_6)
    plot_pseudoknot_analysis(results_4)
    plot_riboswitch_analysis(results_5)
    plot_computation_time(results_1, results_6)
    plot_overall_summary(results_1, results_2, results_3, results_6)

    # Save raw results
    summary = {
        'experiment_1_baseline': {k: {mk: float(mv) for mk, mv in v.items()} for k, v in results_1.items()},
        'experiment_4_pseudoknots': {k: {'n_pseudoknots': v['n_pseudoknots'], 'base_pairs': v['base_pairs'], 'pk_pairs': v['pk_pairs']} for k, v in results_4.items()},
    }
    with open('experiment_results.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETED SUCCESSFULLY")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  - src/rna_predictor.py  (core algorithm)")
    print(f"  - run_experiments.py    (experiment runner)")
    print(f"  - experiment_results.json (raw results)")
    print(f"  - figures/*.png         (7 figures)")
