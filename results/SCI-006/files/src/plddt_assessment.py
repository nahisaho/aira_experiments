"""
Module 1: AlphaFold2 pLDDT-based Docking Suitability Assessment

Evaluates AlphaFold2 predicted structures for molecular docking suitability
based on per-residue pLDDT confidence scores in the binding site region.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class BindingSiteResidue:
    residue_id: int
    residue_name: str
    plddt_score: float
    is_binding_site: bool
    sasa: float  # solvent accessible surface area


class PLDDTAssessor:
    """Assesses docking suitability based on AlphaFold2 pLDDT scores."""

    QUALITY_THRESHOLDS = {
        'very_high': 90.0,
        'confident': 70.0,
        'low': 50.0,
    }

    def __init__(self, plddt_scores: np.ndarray, binding_site_indices: List[int]):
        self.plddt_scores = plddt_scores
        self.binding_site_indices = binding_site_indices
        self.binding_site_plddt = plddt_scores[binding_site_indices]

    def compute_suitability_score(self) -> float:
        """Compute weighted docking suitability score (0-1)."""
        mean_plddt = np.mean(self.binding_site_plddt)
        fraction_confident = np.mean(self.binding_site_plddt > self.QUALITY_THRESHOLDS['confident'])
        min_plddt = np.min(self.binding_site_plddt)

        score = (
            0.4 * (mean_plddt / 100.0) +
            0.4 * fraction_confident +
            0.2 * (min_plddt / 100.0)
        )
        return float(np.clip(score, 0, 1))

    def classify_quality(self) -> str:
        score = self.compute_suitability_score()
        if score >= 0.85:
            return "Excellent - suitable for high-confidence docking"
        elif score >= 0.70:
            return "Good - suitable for docking with minor caveats"
        elif score >= 0.50:
            return "Moderate - docking results should be interpreted cautiously"
        else:
            return "Poor - structure refinement recommended before docking"

    def identify_problematic_residues(self) -> List[int]:
        return [idx for idx in self.binding_site_indices
                if self.plddt_scores[idx] < self.QUALITY_THRESHOLDS['low']]

    def get_statistics(self) -> dict:
        return {
            'mean_plddt_overall': float(np.mean(self.plddt_scores)),
            'mean_plddt_binding_site': float(np.mean(self.binding_site_plddt)),
            'std_plddt_binding_site': float(np.std(self.binding_site_plddt)),
            'min_plddt_binding_site': float(np.min(self.binding_site_plddt)),
            'max_plddt_binding_site': float(np.max(self.binding_site_plddt)),
            'fraction_confident': float(np.mean(self.binding_site_plddt > 70)),
            'fraction_very_high': float(np.mean(self.binding_site_plddt > 90)),
            'suitability_score': self.compute_suitability_score(),
            'quality_class': self.classify_quality(),
            'n_problematic_residues': len(self.identify_problematic_residues()),
        }


def simulate_alphafold_structure(n_residues: int = 300, seed: int = 42) -> Tuple[np.ndarray, List[int]]:
    """Generate simulated AlphaFold2 pLDDT scores for demonstration."""
    rng = np.random.RandomState(seed)

    plddt = np.zeros(n_residues)

    # Divide into structured and loop regions proportionally
    segment_len = n_residues // 6
    for i in range(n_residues):
        region = i // segment_len
        if region % 2 == 0:  # structured
            plddt[i] = rng.normal(88, 5)
        else:  # loop
            plddt[i] = rng.normal(62, 12)

    plddt = np.clip(plddt, 0, 100)

    # Define binding site as ~15% of residues in the middle
    bs_start = n_residues // 4
    bs_size = max(15, n_residues // 7)
    binding_site = list(range(bs_start, min(bs_start + bs_size, n_residues)))

    return plddt, binding_site


def run_plddt_analysis(output_dir: str = "figures"):
    """Run the complete pLDDT analysis and generate figures."""
    print("=" * 60)
    print("Module 1: pLDDT-based Docking Suitability Assessment")
    print("=" * 60)

    # Simulate multiple protein targets
    targets = {
        'CDK2 (kinase)': (300, 42),
        'BRD4 (bromodomain)': (250, 123),
        'SARS-CoV-2 Mpro': (306, 456),
        'PDE5 (phosphodiesterase)': (350, 789),
        'EGFR (kinase)': (280, 321),
    }

    results = {}
    all_plddt_data = {}

    for target_name, (n_res, seed) in targets.items():
        plddt, bs_indices = simulate_alphafold_structure(n_res, seed)
        assessor = PLDDTAssessor(plddt, bs_indices)
        stats = assessor.get_statistics()
        results[target_name] = stats
        all_plddt_data[target_name] = (plddt, bs_indices)

        print(f"\n{target_name}:")
        print(f"  Overall pLDDT: {stats['mean_plddt_overall']:.1f}")
        print(f"  Binding site pLDDT: {stats['mean_plddt_binding_site']:.1f} ± {stats['std_plddt_binding_site']:.1f}")
        print(f"  Suitability Score: {stats['suitability_score']:.3f}")
        print(f"  Quality: {stats['quality_class']}")

    # Figure 1: pLDDT profiles
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    for idx, (name, (plddt, bs_idx)) in enumerate(all_plddt_data.items()):
        ax = axes[idx]
        residues = np.arange(len(plddt))
        ax.plot(residues, plddt, 'b-', alpha=0.5, linewidth=0.8, label='All residues')
        ax.scatter(bs_idx, plddt[bs_idx], c='red', s=8, zorder=5, label='Binding site')
        ax.axhline(y=70, color='orange', linestyle='--', alpha=0.5, label='Confidence threshold')
        ax.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='High confidence')
        ax.set_title(name, fontsize=10)
        ax.set_xlabel('Residue Index')
        ax.set_ylabel('pLDDT')
        ax.set_ylim(0, 105)
        if idx == 0:
            ax.legend(fontsize=7)

    axes[-1].axis('off')
    plt.suptitle('AlphaFold2 pLDDT Profiles with Binding Site Residues', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/plddt_profiles.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Figure 2: Suitability scores comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    names = list(results.keys())
    scores = [results[n]['suitability_score'] for n in names]
    colors = ['green' if s > 0.85 else 'orange' if s > 0.70 else 'red' for s in scores]

    bars = ax1.barh(names, scores, color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Docking Suitability Score')
    ax1.set_title('Docking Suitability Assessment')
    ax1.set_xlim(0, 1)
    ax1.axvline(x=0.85, color='green', linestyle='--', alpha=0.5, label='Excellent')
    ax1.axvline(x=0.70, color='orange', linestyle='--', alpha=0.5, label='Good')
    ax1.legend()
    for bar, score in zip(bars, scores):
        ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                 f'{score:.3f}', va='center', fontsize=9)

    # pLDDT distribution comparison
    bs_means = [results[n]['mean_plddt_binding_site'] for n in names]
    overall_means = [results[n]['mean_plddt_overall'] for n in names]
    x = np.arange(len(names))
    width = 0.35
    ax2.bar(x - width / 2, overall_means, width, label='Overall', color='steelblue', edgecolor='black', linewidth=0.5)
    ax2.bar(x + width / 2, bs_means, width, label='Binding Site', color='coral', edgecolor='black', linewidth=0.5)
    ax2.set_ylabel('Mean pLDDT')
    ax2.set_title('pLDDT: Overall vs Binding Site')
    ax2.set_xticks(x)
    ax2.set_xticklabels([n.split('(')[0].strip() for n in names], rotation=30, ha='right')
    ax2.legend()
    ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/plddt_suitability.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nFigures saved to {output_dir}/")
    return results


if __name__ == '__main__':
    run_plddt_analysis()
