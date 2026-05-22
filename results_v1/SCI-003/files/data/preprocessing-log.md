# Preprocessing Log

- Timestamp: 2026-05-22T12:48:55.624474
- RNA: QC filters (min_genes=200, max_genes=5000, max_mito_pct=20), normalize_total, log1p, HVG, PCA, UMAP
- ATAC: TF-IDF normalization, LSI (TruncatedSVD), UMAP
- Methylation: low-variance CpG check, standardization, PCA, UMAP
- RNA cells retained: 500, genes: 2000
- ATAC cells retained: 500, peaks: 5000
- Methylation cells retained: 500, CpGs: 1000
