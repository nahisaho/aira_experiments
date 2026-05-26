import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('Agg')

import json
import os
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent
SRC_DIR = ROOT / 'src'
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data_utils import (
    AMINO_ACIDS,
    amino_acid_frequency,
    compute_sequence_properties,
    generate_synthetic_cdrh3_dataset,
    pairwise_similarity_matrix,
    summarize_property,
)
from diffusion_model import generate_sequences, save_checkpoint, train_diffusion_model
from optimizer import compute_pareto_front, guided_sampling, select_top_candidates
from pdl1_case_study import compare_with_reference, generate_pdl1_candidates, run_insilico_validation
from property_predictor import evaluate_all_properties, train_property_predictor


def safe_make_figure(filename: Path, title: str, plot_fn: Callable[[], None]) -> None:
    try:
        plt.close('all')
        plot_fn()
        plt.tight_layout()
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        print(f'[FIGURE] Saved {title}: {filename}')
    except Exception as exc:
        print(f'[FIGURE] Failed to generate {title}: {exc}')
    finally:
        plt.close('all')


def sequence_dict_to_lists(items: List[Dict[str, float]], keys: List[str]) -> Dict[str, List[float]]:
    return {key: [float(item[key]) for item in items] for key in keys}


def main() -> None:
    start_time = time.time()
    np.random.seed(42)
    torch.manual_seed(42)
    torch.set_num_threads(max(1, min(4, os.cpu_count() or 1)))

    figures_dir = ROOT / 'figures'
    results_dir = ROOT / 'results'
    figures_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    print('[SETUP] Starting de novo therapeutic antibody design experiment.')
    print('[SETUP] Random seeds set to 42 for NumPy and PyTorch.')

    print('[1/7] Preparing synthetic CDR-H3 datasets...')
    dataset = generate_synthetic_cdrh3_dataset(num_train=1000, num_test=200, seed=42)
    train_sequences = dataset['train_sequences']
    test_sequences = dataset['test_sequences']
    train_properties = dataset['train_properties']
    test_properties = dataset['test_properties']
    print(f'[DATA] Training sequences: {len(train_sequences)} | Test sequences: {len(test_sequences)}')
    print('[DATA] Sequence properties computed for all samples.')

    print('[2/7] Training diffusion model...')
    diffusion_artifacts = train_diffusion_model(
        train_sequences,
        epochs=50,
        batch_size=64,
        lr=1e-4,
        device='cpu',
    )
    for epoch, loss in enumerate(diffusion_artifacts.losses, start=1):
        print(f'[DIFFUSION] Epoch {epoch:02d}/50 - loss={loss:.4f}')
    checkpoint_path = results_dir / 'diffusion_model.pt'
    save_checkpoint(diffusion_artifacts, checkpoint_path)
    print(f'[DIFFUSION] Saved checkpoint to {checkpoint_path}')

    print('[3/7] Training property predictors...')
    property_names = ['binding_affinity', 'stability', 'humanization', 'developability']
    predictors = {}
    for name in property_names:
        print(f'[PREDICTOR] Training {name} predictor...')
        targets = [float(item[name]) for item in train_properties]
        artifacts = train_property_predictor(
            train_sequences,
            targets,
            property_name=name,
            epochs=30,
            batch_size=64,
            lr=1e-3,
            device='cpu',
        )
        predictors[name] = artifacts
        print(f'[PREDICTOR] Finished {name} predictor. Final loss={artifacts.losses[-1]:.4f}')

    print('[4/7] Generating novel CDR-H3 sequences...')
    generated_sequences = generate_sequences(diffusion_artifacts, num_sequences=500, num_steps=8, seed=42, device='cpu')
    generated_predictions = evaluate_all_properties(predictors, generated_sequences, device='cpu')
    generated_physics = compute_sequence_properties(generated_sequences, seed=84)
    for predicted, physical in zip(generated_predictions, generated_physics):
        predicted['expression'] = physical['expression']
        predicted['aggregation'] = physical['aggregation']
        predicted['length'] = physical['length']
    print(f'[GENERATION] Generated {len(generated_sequences)} sequences and evaluated all properties.')

    print('[5/7] Running multi-objective optimization...')
    optimized_candidates, optimization_trajectory = guided_sampling(
        diffusion_artifacts,
        predictors,
        num_candidates=250,
        iterations=8,
        batch_size=50,
        seed=42,
        device='cpu',
    )
    pareto_front = compute_pareto_front(optimized_candidates)
    top_candidates = select_top_candidates(pareto_front, top_k=5)
    top_candidate_sequences = [item['sequence'] for item in top_candidates]
    top_candidate_physics = compute_sequence_properties(top_candidate_sequences, seed=91)
    for item, physical in zip(top_candidates, top_candidate_physics):
        item['expression'] = physical['expression']
        item['aggregation'] = physical['aggregation']
    print(f'[OPTIMIZATION] Optimized pool size: {len(optimized_candidates)} | Pareto front size: {len(pareto_front)}')

    print('[6/7] Running PD-L1 case study...')
    pdl1_candidates = generate_pdl1_candidates(diffusion_artifacts, predictors, num_candidates=40, seed=42, device='cpu')
    validated_pdl1_candidates = run_insilico_validation(pdl1_candidates)
    pdl1_comparison = compare_with_reference(validated_pdl1_candidates)
    print(
        '[PD-L1] Top validated candidate score='
        f"{validated_pdl1_candidates[0]['validated_pdl1_binding']:.3f} | "
        f"Best reference={max(item['pdl1_binding'] for item in pdl1_comparison['references']):.3f}"
    )

    print('[7/7] Generating figures...')
    train_freq = amino_acid_frequency(train_sequences)
    generated_freq = amino_acid_frequency(generated_sequences)
    generated_lengths = [len(seq) for seq in generated_sequences]
    property_arrays = sequence_dict_to_lists(generated_predictions, property_names)

    safe_make_figure(figures_dir / 'training_loss.png', 'training loss curve', lambda: (
        plt.figure(figsize=(7, 4)),
        plt.plot(range(1, len(diffusion_artifacts.losses) + 1), diffusion_artifacts.losses, color='navy', linewidth=2),
        plt.xlabel('Epoch'),
        plt.ylabel('Loss'),
        plt.title('Diffusion Model Training Loss')
    ))

    safe_make_figure(figures_dir / 'generated_length_distribution.png', 'generated length distribution', lambda: (
        plt.figure(figsize=(7, 4)),
        plt.hist(generated_lengths, bins=range(min(generated_lengths), max(generated_lengths) + 2), color='teal', edgecolor='black'),
        plt.xlabel('CDR-H3 Length'),
        plt.ylabel('Count'),
        plt.title('Generated CDR-H3 Length Distribution')
    ))

    def plot_amino_acid_frequency() -> None:
        x = np.arange(len(AMINO_ACIDS))
        width = 0.38
        plt.figure(figsize=(10, 4))
        plt.bar(x - width / 2, [train_freq[aa] for aa in AMINO_ACIDS], width=width, label='Training', color='slateblue')
        plt.bar(x + width / 2, [generated_freq[aa] for aa in AMINO_ACIDS], width=width, label='Generated', color='darkorange')
        plt.xticks(x, AMINO_ACIDS)
        plt.ylabel('Frequency')
        plt.title('Amino Acid Frequency Comparison')
        plt.legend()

    safe_make_figure(figures_dir / 'amino_acid_frequency.png', 'amino acid frequency comparison', plot_amino_acid_frequency)

    def plot_property_distributions() -> None:
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        for ax, name in zip(axes.flatten(), property_names):
            ax.hist(property_arrays[name], bins=20, color='cornflowerblue', edgecolor='black', alpha=0.85)
            ax.set_title(name.replace('_', ' ').title())
            ax.set_xlabel('Score')
            ax.set_ylabel('Count')
        fig.suptitle('Generated Sequence Property Distributions', y=1.02)

    safe_make_figure(figures_dir / 'property_distributions.png', 'property distributions', plot_property_distributions)

    def plot_pareto_front() -> None:
        plt.figure(figsize=(7, 5))
        plt.scatter(
            [item['binding_affinity'] for item in optimized_candidates],
            [item['stability'] for item in optimized_candidates],
            c=[item['humanization'] for item in optimized_candidates],
            cmap='viridis',
            alpha=0.35,
            label='Optimized pool',
        )
        plt.scatter(
            [item['binding_affinity'] for item in pareto_front],
            [item['stability'] for item in pareto_front],
            c=[item['humanization'] for item in pareto_front],
            cmap='viridis',
            edgecolor='black',
            s=60,
            label='Pareto front',
        )
        plt.xlabel('Binding Affinity')
        plt.ylabel('Stability')
        plt.title('Pareto Front of Optimized Candidates')
        plt.colorbar(label='Humanization Score')
        plt.legend()

    safe_make_figure(figures_dir / 'pareto_front.png', 'pareto front', plot_pareto_front)

    def plot_pdl1_binding_scores() -> None:
        top_generated = validated_pdl1_candidates[:5]
        labels = [f'Gen {idx + 1}' for idx in range(len(top_generated))] + [item['name'] for item in pdl1_comparison['references']]
        values = [item['validated_pdl1_binding'] for item in top_generated] + [item['pdl1_binding'] for item in pdl1_comparison['references']]
        colors = ['seagreen'] * len(top_generated) + ['gray'] * len(pdl1_comparison['references'])
        plt.figure(figsize=(9, 4))
        plt.bar(labels, values, color=colors)
        plt.ylabel('PD-L1 Binding Score')
        plt.title('PD-L1 Candidate Scores vs Reference Antibodies')
        plt.xticks(rotation=25, ha='right')

    safe_make_figure(figures_dir / 'pdl1_binding_scores.png', 'PD-L1 binding scores', plot_pdl1_binding_scores)

    def plot_radar_chart() -> None:
        labels = ['binding_affinity', 'stability', 'humanization', 'developability']
        angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
        angles += angles[:1]
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, polar=True)
        for idx, item in enumerate(top_candidates[:5]):
            values = [item[label] for label in labels]
            values += values[:1]
            ax.plot(angles, values, linewidth=2, label=f'C{idx + 1}')
            ax.fill(angles, values, alpha=0.08)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels([label.replace('_', ' ').title() for label in labels])
        ax.set_title('Top Candidate Multi-Property Comparison')
        ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1))

    safe_make_figure(figures_dir / 'multi_property_radar.png', 'multi-property radar chart', plot_radar_chart)

    def plot_sequence_diversity() -> None:
        similarity = pairwise_similarity_matrix(top_candidate_sequences)
        plt.figure(figsize=(6, 5))
        plt.imshow(similarity, cmap='magma', vmin=0, vmax=1)
        plt.colorbar(label='Pairwise Similarity')
        plt.xticks(range(len(top_candidate_sequences)), [f'C{i + 1}' for i in range(len(top_candidate_sequences))])
        plt.yticks(range(len(top_candidate_sequences)), [f'C{i + 1}' for i in range(len(top_candidate_sequences))])
        plt.title('Top Candidate Sequence Diversity')

    safe_make_figure(figures_dir / 'sequence_diversity.png', 'sequence diversity heatmap', plot_sequence_diversity)

    def plot_optimization_trajectory() -> None:
        iterations = [item['iteration'] for item in optimization_trajectory]
        best_scores = [item['best_score'] for item in optimization_trajectory]
        mean_scores = [item['mean_score'] for item in optimization_trajectory]
        plt.figure(figsize=(7, 4))
        plt.plot(iterations, best_scores, marker='o', label='Best score', color='crimson')
        plt.plot(iterations, mean_scores, marker='s', label='Mean score', color='steelblue')
        plt.xlabel('Iteration')
        plt.ylabel('Composite Score')
        plt.title('Optimization Trajectory')
        plt.legend()

    safe_make_figure(figures_dir / 'optimization_trajectory.png', 'optimization trajectory', plot_optimization_trajectory)

    def plot_developability_assessment() -> None:
        plt.figure(figsize=(7, 5))
        plt.scatter(
            [item['expression'] for item in generated_predictions],
            [item['aggregation'] for item in generated_predictions],
            c=[item['developability'] for item in generated_predictions],
            cmap='coolwarm_r',
            alpha=0.65,
        )
        plt.xlabel('Expression')
        plt.ylabel('Aggregation')
        plt.title('Developability Assessment')
        plt.colorbar(label='Developability Score')

    safe_make_figure(figures_dir / 'developability_assessment.png', 'developability assessment', plot_developability_assessment)

    binding_avg, binding_best = summarize_property('binding_affinity', generated_predictions)
    stability_avg, stability_best = summarize_property('stability', generated_predictions)
    humanization_avg, _ = summarize_property('humanization', generated_predictions)

    results = {
        'configuration': {
            'seed': 42,
            'train_sequences': 1000,
            'test_sequences': 200,
            'generated_sequences': 500,
            'diffusion_epochs': 50,
            'predictor_epochs': 30,
        },
        'diffusion_training_loss': diffusion_artifacts.losses,
        'predictor_training_loss': {name: artifacts.losses for name, artifacts in predictors.items()},
        'summary': {
            'num_generated_sequences': len(generated_sequences),
            'average_binding_affinity': binding_avg,
            'best_binding_affinity': binding_best,
            'average_stability': stability_avg,
            'best_stability': stability_best,
            'average_humanization_score': humanization_avg,
            'num_pareto_optimal_candidates': len(pareto_front),
            'runtime_seconds': time.time() - start_time,
        },
        'top_candidates': top_candidates,
        'pareto_front': pareto_front[:25],
        'optimization_trajectory': optimization_trajectory,
        'pdl1_case_study': {
            'top_validated_candidates': validated_pdl1_candidates[:10],
            'reference_comparison': pdl1_comparison,
        },
    }

    results_path = results_dir / 'experiment_results.json'
    with results_path.open('w', encoding='utf-8') as handle:
        json.dump(results, handle, indent=2)
    print(f'[RESULTS] Saved experiment results to {results_path}')

    print('[SUMMARY] Number of generated sequences:', len(generated_sequences))
    print(f'[SUMMARY] Average / best binding affinity: {binding_avg:.3f} / {binding_best:.3f}')
    print(f'[SUMMARY] Average / best stability: {stability_avg:.3f} / {stability_best:.3f}')
    print(f'[SUMMARY] Average humanization score: {humanization_avg:.3f}')
    print(f'[SUMMARY] Number of Pareto-optimal candidates: {len(pareto_front)}')
    print(
        '[SUMMARY] PD-L1 case study best candidate vs best reference: '
        f"{validated_pdl1_candidates[0]['validated_pdl1_binding']:.3f} vs "
        f"{max(item['pdl1_binding'] for item in pdl1_comparison['references']):.3f}"
    )
    print(f'[DONE] Total runtime: {time.time() - start_time:.2f} seconds')


if __name__ == '__main__':
    main()
