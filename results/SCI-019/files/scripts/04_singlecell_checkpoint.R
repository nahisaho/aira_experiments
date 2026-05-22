#!/usr/bin/env Rscript

set.seed(42)

run_singlecell_checkpoint <- function(base_dir = NULL) {
  message("[04] Starting single-cell checkpoint analysis...")

  suppressPackageStartupMessages({
    pkgs <- c("Seurat", "ggplot2", "dplyr", "patchwork", "viridis", "RColorBrewer", "scales")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  })

  resolve_base_dir <- function(path_hint = NULL) {
    if (!is.null(path_hint)) return(normalizePath(path_hint, mustWork = TRUE))
    wd <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(wd) == "scripts") dirname(wd) else wd
  }

  base_dir <- resolve_base_dir(base_dir)
  fig_dir <- file.path(base_dir, "figures")
  res_dir <- file.path(base_dir, "results")
  data_dir <- file.path(base_dir, "data")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  checkpoint_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "PDCD1LG2", "IDO1", "FOXP3")
  marker_genes <- c(
    "IL7R", "LTB", "MALAT1", "TRAC", "CD3D", "NKG7", "KLRD1", "GNLY",
    "MS4A1", "CD79A", "HLA-DRA", "FCER1A", "CLEC10A", "LST1", "S100A8", "TYROBP"
  )
  other_genes <- sprintf("GENE_%04d", seq_len(2000 - length(unique(c(checkpoint_genes, marker_genes)))))
  genes <- unique(c(checkpoint_genes, marker_genes, other_genes))
  genes <- genes[seq_len(2000)]

  cell_types <- c("CD4 T", "CD8 T", "B cells", "NK cells", "Monocytes", "DCs")
  cell_type_prob <- c(0.30, 0.25, 0.15, 0.10, 0.12, 0.08)
  conditions <- rep(c("RA_active", "HC"), each = 2500)
  cell_ids <- sprintf("Cell_%04d", seq_len(5000))
  cell_type_assign <- sample(cell_types, 5000, replace = TRUE, prob = cell_type_prob)
  meta <- data.frame(
    cell_id = cell_ids,
    condition = conditions,
    cell_type = cell_type_assign,
    stringsAsFactors = FALSE
  )
  meta$celltype_condition <- paste(meta$cell_type, meta$condition, sep = " | ")
  rownames(meta) <- meta$cell_id

  base_mean <- rgamma(length(genes), shape = 2, rate = 1)
  names(base_mean) <- genes
  celltype_effects <- matrix(0, nrow = length(genes), ncol = length(cell_types), dimnames = list(genes, cell_types))
  condition_effects <- matrix(0, nrow = length(genes), ncol = 2, dimnames = list(genes, c("RA_active", "HC")))

  celltype_effects[c("IL7R", "LTB", "TRAC", "CD3D"), "CD4 T"] <- 2.3
  celltype_effects[c("TRAC", "CD3D", "NKG7"), "CD8 T"] <- 2.1
  celltype_effects[c("MS4A1", "CD79A"), "B cells"] <- 2.5
  celltype_effects[c("NKG7", "KLRD1", "GNLY"), "NK cells"] <- 2.8
  celltype_effects[c("LST1", "S100A8", "TYROBP"), "Monocytes"] <- 2.6
  celltype_effects[c("HLA-DRA", "FCER1A", "CLEC10A"), "DCs"] <- 2.7

  celltype_effects[c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"), c("CD4 T", "CD8 T")] <- 1.4
  celltype_effects[c("CD274", "PDCD1LG2", "IDO1"), c("Monocytes", "DCs", "B cells")] <- 1.5
  celltype_effects["FOXP3", "CD4 T"] <- 1.9

  condition_effects[c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "FOXP3"), "RA_active"] <- 0.85
  condition_effects[c("CD274", "PDCD1LG2", "IDO1"), "RA_active"] <- 0.95
  condition_effects[c("IL7R", "TRAC"), "HC"] <- 0.25

  exhausted_signal <- rep(0, 5000)
  exhausted_signal[meta$condition == "RA_active" & meta$cell_type %in% c("CD4 T", "CD8 T")] <- rnorm(sum(meta$condition == "RA_active" & meta$cell_type %in% c("CD4 T", "CD8 T")), mean = 0.9, sd = 0.2)
  exhausted_signal[meta$condition == "HC" & meta$cell_type %in% c("CD4 T", "CD8 T")] <- rnorm(sum(meta$condition == "HC" & meta$cell_type %in% c("CD4 T", "CD8 T")), mean = 0.2, sd = 0.12)

  count_matrix <- matrix(0L, nrow = length(genes), ncol = 5000, dimnames = list(genes, cell_ids))
  checkpoint_index <- match(checkpoint_genes, genes)

  for (i in seq_len(5000)) {
    ct <- meta$cell_type[i]
    cond <- meta$condition[i]
    mu <- base_mean * exp(celltype_effects[, ct] + condition_effects[, cond] + rnorm(length(genes), 0, 0.15))
    if (ct %in% c("CD4 T", "CD8 T")) {
      mu[checkpoint_index[c(1, 4, 5, 6)]] <- mu[checkpoint_index[c(1, 4, 5, 6)]] * exp(exhausted_signal[i])
    }
    if (ct %in% c("Monocytes", "DCs", "B cells") && cond == "RA_active") {
      mu[checkpoint_index[c(2, 7, 8)]] <- mu[checkpoint_index[c(2, 7, 8)]] * 1.8
    }
    if (ct == "CD4 T" && cond == "RA_active") {
      mu[checkpoint_index[9]] <- mu[checkpoint_index[9]] * 1.35
    }
    count_matrix[, i] <- rnbinom(length(mu), mu = pmax(mu, 0.01), size = 1.2)
  }

  saveRDS(list(counts = count_matrix, metadata = meta), file.path(data_dir, "synthetic_scrnaseq_checkpoint.rds"))

  seu <- Seurat::CreateSeuratObject(counts = count_matrix, meta.data = meta, project = "AutoimmuneCheckpoint")
  seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- Seurat::RunUMAP(seu, dims = 1:20, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:20, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = 0.55, verbose = FALSE)

  cluster_summary <- seu@meta.data %>%
    dplyr::count(seurat_clusters, cell_type, condition, name = "n_cells") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(cluster_fraction = n_cells / sum(n_cells)) %>%
    dplyr::ungroup()
  write.csv(cluster_summary, file.path(res_dir, "checkpoint_cluster_summary.csv"), row.names = FALSE)

  p_celltype <- Seurat::DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE) +
    ggplot2::labs(title = "UMAP by cell type")
  p_condition <- Seurat::DimPlot(seu, reduction = "umap", group.by = "condition") +
    ggplot2::labs(title = "UMAP by condition")
  ggplot2::ggsave(file.path(fig_dir, "scrnaseq_umap_overview.pdf"), p_celltype + p_condition, width = 12, height = 5)

  feature_panel <- Seurat::FeaturePlot(seu, reduction = "umap", features = c("PDCD1", "HAVCR2", "CD274", "FOXP3"), ncol = 2) &
    ggplot2::theme(legend.position = "right")
  ggplot2::ggsave(file.path(fig_dir, "scrnaseq_checkpoint_featureplots.pdf"), feature_panel, width = 12, height = 9)

  violin_panel <- Seurat::VlnPlot(seu, features = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "CD274", "IDO1"), group.by = "cell_type", pt.size = 0, ncol = 3) &
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(fig_dir, "scrnaseq_checkpoint_violins.pdf"), violin_panel, width = 14, height = 8)

  dotplot <- Seurat::DotPlot(seu, features = checkpoint_genes, group.by = "celltype_condition") +
    ggplot2::scale_color_viridis_c(option = "viridis") +
    ggplot2::labs(title = "Checkpoint expression by cell type and condition") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(fig_dir, "scrnaseq_checkpoint_dotplot.pdf"), dotplot, width = 13, height = 6)

  expr_mat <- Seurat::GetAssayData(seu, slot = "data")
  tcell_idx <- seu$cell_type %in% c("CD4 T", "CD8 T")
  exhausted_flag <- expr_mat["PDCD1", ] > 1.2 & expr_mat["HAVCR2", ] > 1.0 & tcell_idx
  exhausted_freq <- seu@meta.data %>%
    dplyr::mutate(exhausted = exhausted_flag[rownames(seu@meta.data)]) %>%
    dplyr::filter(cell_type %in% c("CD4 T", "CD8 T")) %>%
    dplyr::group_by(condition, cell_type) %>%
    dplyr::summarise(exhausted_frequency = mean(exhausted), n_cells = dplyr::n(), .groups = "drop")
  write.csv(exhausted_freq, file.path(res_dir, "checkpoint_exhausted_tcell_frequency.csv"), row.names = FALSE)

  de_results <- lapply(unique(seu$cell_type), function(ct) {
    subset_obj <- subset(seu, subset = cell_type == ct)
    Idents(subset_obj) <- subset_obj$condition
    mk <- Seurat::FindMarkers(
      subset_obj,
      ident.1 = "RA_active",
      ident.2 = "HC",
      features = checkpoint_genes,
      logfc.threshold = 0,
      min.pct = 0.05,
      test.use = "wilcox",
      verbose = FALSE
    )
    mk$gene <- rownames(mk)
    mk$cell_type <- ct
    mk
  })
  de_df <- dplyr::bind_rows(de_results)
  de_df$padj <- p.adjust(de_df$p_val, method = "BH")
  write.csv(de_df, file.path(res_dir, "checkpoint_DE.csv"), row.names = FALSE)

  checkpoint_summary <- seu@meta.data %>%
    dplyr::mutate(
      PDCD1 = as.numeric(expr_mat["PDCD1", rownames(seu@meta.data)]),
      HAVCR2 = as.numeric(expr_mat["HAVCR2", rownames(seu@meta.data)]),
      CD274 = as.numeric(expr_mat["CD274", rownames(seu@meta.data)]),
      FOXP3 = as.numeric(expr_mat["FOXP3", rownames(seu@meta.data)])
    ) %>%
    dplyr::group_by(condition, cell_type) %>%
    dplyr::summarise(
      PDCD1_mean = mean(PDCD1),
      HAVCR2_mean = mean(HAVCR2),
      CD274_mean = mean(CD274),
      FOXP3_mean = mean(FOXP3),
      .groups = "drop"
    )
  write.csv(checkpoint_summary, file.path(res_dir, "checkpoint_expression_summary.csv"), row.names = FALSE)

  message("[04] Single-cell checkpoint analysis completed.")
  invisible(list(seurat = seu, de = de_df, exhausted = exhausted_freq, summary = checkpoint_summary))
}

if (sys.nframe() == 0L) {
  run_singlecell_checkpoint()
}
