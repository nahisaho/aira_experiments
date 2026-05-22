# Preprocessing log

- Seed policy: all analysis scripts initialize `set.seed(42)`.
- Multi-omics transcriptome: negative binomial counts are normalized with edgeR log-CPM before limma analysis.
- Multi-omics proteome/metabolome: abundance matrices are simulated on continuous scales and exported directly.
- Immune deconvolution: LM22-like proportions are sampled with group-specific concentration shifts and normalized to sum to 1 per sample.
- Single-cell RNA-seq: raw integer counts are simulated first, then processed through Seurat normalization, scaling, PCA, UMAP, and clustering.
- Treatment prediction: clinical, genomic, immune, and cytokine features are concatenated into a patient-level modeling table.
- Tolerance model: ODE trajectories are summarized into endpoint efficacy scores and Monte Carlo robustness metrics.
