# Experimental Protocol

## Objective

Construct a synthetic but biologically plausible heterogeneous knowledge graph and compare TransE, RotatE, and ComplEx for link prediction and COVID-19 drug repurposing.

## Materials and software

- Python 3
- NumPy, pandas, matplotlib, networkx
- PyTorch if available; NumPy fallback otherwise
- Fixed seeds: random=42, numpy=42, torch=42

## Graph schema

Node types: Drug, Disease, Gene, Pathway, Phenotype.
Relation types: treats, targets, associated_with, part_of, manifests_as, interacts_with.
The graph includes curated COVID-19 entities and synthetic background entities with controlled noise.

## Procedure

1. Generate 200+ entities with typed identifiers and descriptive names.
2. Insert 500+ edges using biologically plausible templates plus random noisy links.
3. Reserve a subset of drug-disease edges for validation/testing, ensuring COVID-19 candidate edges are excluded from training when needed.
4. Train TransE, RotatE, and ComplEx for 60 epochs with negative sampling.
5. Perform 3-fold cross-validation on train/validation/test partitions.
6. Compute MRR, Hits@1, Hits@3, Hits@10, and AUROC for each fold.
7. Compare against DegreePrior and Random baselines.
8. Score all unseen drug-COVID-19 pairs and retain top 10 candidates.
9. Generate shortest-path explanations through gene/pathway intermediates.
10. Plot loss curves, metric comparisons, candidate ranking, and COVID-19 subgraph structure.

## Randomization and controls

- Triples are shuffled before splitting.
- Each fold uses a deterministic seed offset to preserve reproducibility.
- Random baseline uses the same candidate pool.
- DegreePrior baseline ranks drugs by graph proximity and target overlap counts.

## Validation strategy

Internal validation uses 3-fold cross-validation. External validation is theoretical only because the graph is synthetic. To compensate, the case study is anchored around real COVID-19 drugs and mechanistic intermediates. Because the system has multiple components (graph construction, embedding model, explanation layer), an ablation-style comparison is included via: (i) random ranking baseline, (ii) degree-based heuristic baseline, and (iii) each embedding model individually.

## Synthetic data assumptions and limitations

The graph is not a licensed biomedical KG and should not be interpreted as clinical evidence. It encodes realistic motifs, not verified causal truth. Noise is intentionally added to avoid trivial predictions and emulate curation incompleteness.
