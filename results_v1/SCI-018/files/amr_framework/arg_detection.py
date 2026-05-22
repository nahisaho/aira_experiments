from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import DATA_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything


@dataclass
class ARGRecord:
    gene_family: str
    category: str
    sequence: str


class ARGDatabase:
    """Synthetic CARD-like ARG database with representative gene families."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.catalog = [
            ("blaTEM", "beta-lactamase"),
            ("blaCTX-M", "beta-lactamase"),
            ("blaOXA", "beta-lactamase"),
            ("aac(3)-II", "aminoglycoside"),
            ("aadA", "aminoglycoside"),
            ("aph(3')-Ia", "aminoglycoside"),
            ("tetA", "tetracycline"),
            ("tetM", "tetracycline"),
            ("tetX", "tetracycline"),
            ("ermB", "macrolide"),
            ("mphA", "macrolide"),
            ("msrD", "macrolide"),
            ("qnrS", "fluoroquinolone"),
            ("qepA", "fluoroquinolone"),
            ("oqxAB", "fluoroquinolone"),
            ("sul1", "sulfonamide"),
            ("sul2", "sulfonamide"),
            ("dfrA", "trimethoprim"),
            ("vanA", "glycopeptide"),
            ("mcr-1", "polymyxin"),
        ]
        self.records = self._build_database()

    @staticmethod
    def _stable_seed(name: str) -> int:
        return sum((idx + 1) * ord(ch) for idx, ch in enumerate(name))

    def _generate_gene_sequence(self, name: str, category: str, length: int = 420) -> str:
        bases = np.array(list("ACGT"))
        rng = np.random.default_rng(self.seed + self._stable_seed(name + category))
        seq = rng.choice(bases, size=length)
        category_motifs = {
            "beta-lactamase": "ATGGCCTACGTTACCGGATCC",
            "aminoglycoside": "ATGAAGGCGTTCGACCTGACC",
            "tetracycline": "ATGTTCCGCGACTACGATGGC",
            "macrolide": "ATGCGTACCGAGTTCGACTTC",
            "fluoroquinolone": "ATGGACGTTGCCGACTTCGAC",
            "sulfonamide": "ATGCCGGTTACCGACTACGGA",
            "trimethoprim": "ATGACCTTCGACGGTTACTCC",
            "glycopeptide": "ATGGTCGACCTTACCGACTTC",
            "polymyxin": "ATGCGACCTTGGTACCGATTC",
        }
        motif = np.array(list(category_motifs[category]))
        start = int(rng.integers(40, length - len(motif) - 40))
        seq[start : start + len(motif)] = motif
        return "".join(seq.tolist())

    def _build_database(self) -> dict[str, ARGRecord]:
        return {
            gene: ARGRecord(gene, category, self._generate_gene_sequence(gene, category))
            for gene, category in self.catalog
        }

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"gene_family": gene, "category": rec.category, "sequence_length": len(rec.sequence)}
                for gene, rec in self.records.items()
            ]
        )


class KmerAligner:
    def __init__(self, k: int = 21) -> None:
        self.k = k

    def _kmers(self, sequence: str) -> list[str]:
        if len(sequence) < self.k:
            return [sequence]
        return [sequence[idx : idx + self.k] for idx in range(len(sequence) - self.k + 1)]

    def build_index(self, target: str) -> dict[str, list[int]]:
        index: dict[str, list[int]] = defaultdict(list)
        for pos, kmer in enumerate(self._kmers(target)):
            index[kmer].append(pos)
        return index

    def align(self, query: str, target: str, target_index: dict[str, list[int]] | None = None) -> dict[str, Any]:
        query_kmers = self._kmers(query)
        index = target_index or self.build_index(target)
        shared = 0
        offset_votes: Counter[int] = Counter()
        for qpos, kmer in enumerate(query_kmers):
            positions = index.get(kmer, [])
            if positions:
                shared += 1
                for tpos in positions[:4]:
                    offset_votes[tpos - qpos] += 1
        coverage = shared / max(1, len(query_kmers))
        if not offset_votes:
            return {"coverage": coverage, "identity": 0.0, "best_start": None, "shared_kmers": shared}

        best_identity = 0.0
        best_start = None
        for offset, _votes in offset_votes.most_common(6):
            start = max(0, min(len(target) - len(query), offset))
            segment = target[start : start + len(query)]
            if len(segment) != len(query):
                continue
            identity = sum(a == b for a, b in zip(query, segment)) / len(query)
            if identity > best_identity:
                best_identity = identity
                best_start = start
        return {
            "coverage": coverage,
            "identity": best_identity,
            "best_start": best_start,
            "shared_kmers": shared,
        }


class ARGDetectionPipeline:
    def __init__(self, seed: int = 42) -> None:
        seed_everything(seed)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.database = ARGDatabase(seed=seed)
        self.aligner = KmerAligner(k=21)

    def _random_dna(self, length: int) -> str:
        return "".join(self.rng.choice(np.array(list("ACGT")), size=length).tolist())

    def _mutate_sequence(self, sequence: str, mutation_rate: float) -> str:
        bases = np.array(list("ACGT"))
        mutated = list(sequence)
        for idx, base in enumerate(mutated):
            if self.rng.random() < mutation_rate:
                choices = bases[bases != base]
                mutated[idx] = str(self.rng.choice(choices))
        return "".join(mutated)

    def simulate_genome(self, n_genes: int = 6, mutation_rate: float = 0.02) -> tuple[str, list[dict[str, Any]]]:
        n_genes = max(1, min(n_genes, len(self.database.records)))
        selected = self.rng.choice(list(self.database.records.keys()), size=n_genes, replace=False)
        genome_parts = [self._random_dna(500)]
        inserted: list[dict[str, Any]] = []
        for gene_name in selected:
            record = self.database.records[str(gene_name)]
            mutated = self._mutate_sequence(record.sequence, mutation_rate)
            genome_parts.append(mutated)
            genome_parts.append(self._random_dna(int(self.rng.integers(80, 150))))
            inserted.append(
                {
                    "gene_family": record.gene_family,
                    "category": record.category,
                    "mutation_rate": mutation_rate,
                    "reference_length": len(record.sequence),
                }
            )
        genome_parts.append(self._random_dna(600))
        return "".join(genome_parts), inserted

    def detect_args(self, genome_seq: str, genome_id: str = "genome_0") -> pd.DataFrame:
        genome_index = self.aligner.build_index(genome_seq)
        detections = []
        for gene_name, record in self.database.records.items():
            alignment = self.aligner.align(record.sequence, genome_seq, genome_index)
            presence = int(alignment["coverage"] >= 0.60 and alignment["identity"] >= 0.75)
            detections.append(
                {
                    "genome_id": genome_id,
                    "gene_family": gene_name,
                    "category": record.category,
                    "presence": presence,
                    "coverage": round(float(alignment["coverage"]), 4),
                    "identity": round(float(alignment["identity"]), 4),
                    "shared_kmers": int(alignment["shared_kmers"]),
                    "best_start": alignment["best_start"],
                }
            )
        return pd.DataFrame(detections)

    def batch_process(self, n_genomes: int = 10) -> pd.DataFrame:
        results = []
        genome_records = []
        for idx in range(n_genomes):
            n_genes = int(self.rng.integers(4, 9))
            mutation_rate = float(self.rng.uniform(0.005, 0.05))
            genome_seq, inserted = self.simulate_genome(n_genes=n_genes, mutation_rate=mutation_rate)
            genome_id = f"genome_{idx + 1:02d}"
            detections = self.detect_args(genome_seq, genome_id=genome_id)
            results.append(detections)
            genome_records.append(
                {
                    "genome_id": genome_id,
                    "genome_length": len(genome_seq),
                    "n_inserted_args": len(inserted),
                    "inserted_genes": ";".join(rec["gene_family"] for rec in inserted),
                    "mutation_rate": mutation_rate,
                }
            )
        detection_df = pd.concat(results, ignore_index=True)
        detection_path = RESULTS_DIR / "arg_detections.csv"
        detection_df.to_csv(detection_path, index=False)
        genome_df = pd.DataFrame(genome_records)
        genome_path = DATA_DIR / "synthetic_genomes.csv"
        genome_df.to_csv(genome_path, index=False)
        db_path = DATA_DIR / "mock_card_database.csv"
        self.database.to_dataframe().to_csv(db_path, index=False)
        log_event(
            phase="component_1",
            event_type="file_written",
            skill_or_tool="ARGDetectionPipeline",
            files_written=[
                str(detection_path.relative_to(ROOT)),
                str(genome_path.relative_to(ROOT)),
                str(db_path.relative_to(ROOT)),
            ],
            handoff_out={"n_genomes": n_genomes, "detections": int(detection_df["presence"].sum())},
        )
        return detection_df


def run_component(seed: int = 42, n_genomes: int = 12) -> dict[str, Any]:
    pipeline = ARGDetectionPipeline(seed=seed)
    detections = pipeline.batch_process(n_genomes=n_genomes)
    presence_matrix = detections.pivot_table(
        index="genome_id", columns="gene_family", values="presence", fill_value=0
    )
    matrix_path = RESULTS_DIR / "arg_presence_matrix.csv"
    presence_matrix.to_csv(matrix_path)
    summary = {
        "n_genomes": int(detections["genome_id"].nunique()),
        "n_gene_families": int(detections["gene_family"].nunique()),
        "positive_detections": int(detections["presence"].sum()),
        "mean_identity_positive": round(float(detections.loc[detections["presence"] == 1, "identity"].mean()), 4),
        "mean_coverage_positive": round(float(detections.loc[detections["presence"] == 1, "coverage"].mean()), 4),
        "top_categories": detections.groupby("category")["presence"].sum().sort_values(ascending=False).head(5).to_dict(),
    }
    save_json(RESULTS_DIR / "arg_detection_summary.json", summary)
    log_event(
        phase="component_1",
        event_type="handoff_completed",
        skill_or_tool="arg_detection",
        handoff_out=summary,
        files_written=[str(matrix_path.relative_to(ROOT))],
    )
    return summary
