"""
Module 5: Activity Cliff Detection and Chemical Space Exploration

Detects activity cliffs in molecular datasets and implements chemical
space exploration strategies using molecular fingerprints and similarity.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist, squareform
from dataclasses import dataclass
from typing import List, Tuple, Dict


@dataclass
class Molecule:
    mol_id: str
    smiles: str
    activity: float  # pIC50
    fingerprint: np.ndarray
    cluster_id: int = -1
    is_cliff: bool = False


@dataclass
class ActivityCliff:
    mol_a: str
    mol_b: str
    similarity: float
    activity_diff: float
    cliff_score: float


class ActivityCliffDetector:
    """Detect activity cliffs based on structural similarity and activity difference."""

    def __init__(self, similarity_threshold: float = 0.8,
                 activity_threshold: float = 2.0):
        self.sim_threshold = similarity_threshold
        self.act_threshold = activity_threshold

    def compute_tanimoto_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        intersection = np.sum(np.minimum(fp1, fp2))
        union = np.sum(np.maximum(fp1, fp2))
        return intersection / max(union, 1e-10)

    def detect_cliffs(self, molecules: List[Molecule]) -> List[ActivityCliff]:
        cliffs = []
        n = len(molecules)

        for i in range(n):
            for j in range(i + 1, n):
                sim = self.compute_tanimoto_similarity(
                    molecules[i].fingerprint, molecules[j].fingerprint
                )
                act_diff = abs(molecules[i].activity - molecules[j].activity)

                if sim >= self.sim_threshold and act_diff >= self.act_threshold:
                    cliff_score = act_diff * sim
                    cliffs.append(ActivityCliff(
                        mol_a=molecules[i].mol_id,
                        mol_b=molecules[j].mol_id,
                        similarity=sim,
                        activity_diff=act_diff,
                        cliff_score=cliff_score
                    ))
                    molecules[i].is_cliff = True
                    molecules[j].is_cliff = True

        return sorted(cliffs, key=lambda c: c.cliff_score, reverse=True)


class ChemicalSpaceExplorer:
    """Explore chemical space using dimensionality reduction and clustering."""

    def __init__(self, n_clusters: int = 5):
        self.n_clusters = n_clusters

    def reduce_dimensions(self, fingerprints: np.ndarray, perplexity: int = 30) -> np.ndarray:
        tsne = TSNE(n_components=2, perplexity=min(perplexity, len(fingerprints) - 1),
                     random_state=42, n_iter=1000)
        return tsne.fit_transform(fingerprints)

    def cluster_molecules(self, fingerprints: np.ndarray) -> np.ndarray:
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
        return kmeans.fit_predict(fingerprints)

    def compute_diversity(self, fingerprints: np.ndarray) -> float:
        if len(fingerprints) < 2:
            return 0.0
        distances = pdist(fingerprints, metric='jaccard')
        return float(np.mean(distances))


def generate_molecular_dataset(n_molecules: int = 200, fp_size: int = 128,
                                seed: int = 42) -> List[Molecule]:
    """Generate a synthetic molecular dataset with activity cliffs."""
    rng = np.random.RandomState(seed)
    molecules = []

    for i in range(n_molecules):
        # Generate fingerprint
        fp = (rng.random(fp_size) > 0.7).astype(float)

        # Base activity from fingerprint
        activity = 5.0 + np.sum(fp[:20]) * 0.3 + rng.normal(0, 0.5)

        # Introduce activity cliffs for some molecules
        if i > 0 and i % 15 == 0:
            # Similar to previous molecule but different activity
            fp = molecules[i - 1].fingerprint.copy()
            fp[rng.randint(0, fp_size, 3)] = 1 - fp[rng.randint(0, fp_size, 3)]
            activity = molecules[i - 1].activity + rng.choice([-3, 3]) + rng.normal(0, 0.3)

        activity = np.clip(activity, 3.0, 11.0)
        smiles = f"C{'C' * rng.randint(3, 10)}O{'N' * rng.randint(0, 3)}"

        molecules.append(Molecule(
            mol_id=f"MOL-{i + 1:04d}",
            smiles=smiles,
            activity=activity,
            fingerprint=fp
        ))

    return molecules


def run_activity_cliff_analysis(output_dir: str = "figures"):
    """Run activity cliff detection and chemical space exploration."""
    print("=" * 60)
    print("Module 5: Activity Cliff Detection & Chemical Space Exploration")
    print("=" * 60)

    # Generate dataset
    molecules = generate_molecular_dataset(n_molecules=200, seed=42)
    print(f"Dataset: {len(molecules)} molecules")

    # Detect activity cliffs
    detector = ActivityCliffDetector(similarity_threshold=0.75, activity_threshold=1.5)
    cliffs = detector.detect_cliffs(molecules)

    print(f"\nActivity Cliffs Detected: {len(cliffs)}")
    print(f"Molecules involved in cliffs: {sum(1 for m in molecules if m.is_cliff)}")

    if cliffs:
        print(f"\nTop 5 Activity Cliffs:")
        for i, cliff in enumerate(cliffs[:5]):
            print(f"  {i + 1}. {cliff.mol_a} ↔ {cliff.mol_b}: "
                  f"Sim={cliff.similarity:.3f}, ΔpIC50={cliff.activity_diff:.2f}, "
                  f"Score={cliff.cliff_score:.2f}")

    # Chemical space exploration
    explorer = ChemicalSpaceExplorer(n_clusters=6)
    fingerprints = np.array([m.fingerprint for m in molecules])
    coords_2d = explorer.reduce_dimensions(fingerprints)
    cluster_labels = explorer.cluster_molecules(fingerprints)
    diversity = explorer.compute_diversity(fingerprints)

    for mol, cl in zip(molecules, cluster_labels):
        mol.cluster_id = cl

    print(f"\nChemical Space Diversity: {diversity:.3f}")
    print(f"Clusters: {len(set(cluster_labels))}")

    # Figure 6: Activity cliff and chemical space
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Chemical space colored by activity
    ax = axes[0, 0]
    activities = [m.activity for m in molecules]
    sc = ax.scatter(coords_2d[:, 0], coords_2d[:, 1], c=activities,
                    cmap='RdYlGn', s=30, alpha=0.7, edgecolors='gray', linewidths=0.3)
    plt.colorbar(sc, ax=ax, label='pIC50')
    ax.set_title('Chemical Space (colored by pIC50)')
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')

    # Chemical space colored by cluster
    ax = axes[0, 1]
    sc = ax.scatter(coords_2d[:, 0], coords_2d[:, 1], c=cluster_labels,
                    cmap='Set2', s=30, alpha=0.7, edgecolors='gray', linewidths=0.3)
    # Highlight cliff molecules
    cliff_mask = [m.is_cliff for m in molecules]
    if any(cliff_mask):
        ax.scatter(coords_2d[cliff_mask, 0], coords_2d[cliff_mask, 1],
                   facecolors='none', edgecolors='red', s=80, linewidths=2,
                   label='Activity cliff')
        ax.legend()
    ax.set_title('Chemical Space (clusters + cliffs)')
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')

    # Activity distribution
    ax = axes[1, 0]
    cliff_activities = [m.activity for m in molecules if m.is_cliff]
    non_cliff_activities = [m.activity for m in molecules if not m.is_cliff]
    ax.hist(non_cliff_activities, bins=25, alpha=0.6, label='Non-cliff', color='steelblue', density=True)
    if cliff_activities:
        ax.hist(cliff_activities, bins=15, alpha=0.6, label='Cliff', color='red', density=True)
    ax.set_xlabel('pIC50')
    ax.set_ylabel('Density')
    ax.set_title('Activity Distribution')
    ax.legend()

    # SALI (Structure-Activity Landscape Index)
    ax = axes[1, 1]
    if cliffs:
        sims = [c.similarity for c in cliffs]
        diffs = [c.activity_diff for c in cliffs]
        scores = [c.cliff_score for c in cliffs]
        sc = ax.scatter(sims, diffs, c=scores, cmap='hot_r', s=40, alpha=0.7,
                        edgecolors='black', linewidths=0.3)
        plt.colorbar(sc, ax=ax, label='Cliff Score')
    ax.set_xlabel('Tanimoto Similarity')
    ax.set_ylabel('|ΔpIC50|')
    ax.set_title('Structure-Activity Landscape Index (SALI)')
    ax.axhline(y=1.5, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0.75, color='gray', linestyle='--', alpha=0.5)

    plt.suptitle('Activity Cliff Detection & Chemical Space Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/activity_cliffs.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nFigures saved to {output_dir}/")

    return {
        'n_cliffs': len(cliffs),
        'n_cliff_molecules': sum(1 for m in molecules if m.is_cliff),
        'diversity': diversity,
        'n_clusters': len(set(cluster_labels)),
        'top_cliff_score': cliffs[0].cliff_score if cliffs else 0.0,
    }


if __name__ == '__main__':
    run_activity_cliff_analysis()
