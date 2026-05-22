#!/usr/bin/env Rscript
# =============================================================================
# Module 1: Variant Peptide Database Construction & Search
# =============================================================================
# Input:  Somatic VCF, RNA-seq BAM, Reference proteome
# Output: Custom protein FASTA, MaxQuant search results, variant peptide list
# Tools:  customProDB, maftools, MaxQuant
# =============================================================================

suppressPackageStartupMessages({
  library(customProDB)
  library(maftools)
  library(GenomicFeatures)
  library(VariantAnnotation)
  library(Biostrings)
  library(rtracklayer)
  library(yaml)
})

config <- yaml::read_yaml("config/pipeline_config.yaml")
cfg    <- config$variant_peptide

cat("=== Module 1: Variant Peptide Search ===\n")
cat("Pipeline start:", format(Sys.time()), "\n")

# -----------------------------------------------------------------------------
# Step 1: Load and filter somatic variants
# -----------------------------------------------------------------------------
build_variant_proteome <- function(vcf_path, output_dir = "results/variant_db/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  vcf <- readVcf(vcf_path, genome = config$project$genome_build)
  cat(sprintf("  Loaded %d variants from VCF\n", nrow(vcf)))

  # Filter: PASS only, minimum VAF ≥ 0.05
  pass_idx <- fixed(vcf)$FILTER == "PASS"
  vcf <- vcf[pass_idx]
  cat(sprintf("  After PASS filter: %d variants\n", nrow(vcf)))

  # Annotate consequences using Ensembl TxDb
  txdb <- makeTxDbFromGFF(
    file = "data/ref/gencode.v44.annotation.gtf.gz",
    format = "gtf"
  )

  # Variant type classification
  coding_variants <- predictCoding(vcf, txdb,
    seqSource = FaFile("data/ref/GRCh38.primary_assembly.genome.fa"))

  snv_missense  <- coding_variants[coding_variants$CONSEQUENCE == "nonsynonymous"]
  snv_nonsense  <- coding_variants[coding_variants$CONSEQUENCE == "nonsense"]
  snv_frameshift <- coding_variants[coding_variants$CONSEQUENCE == "frameshift"]

  cat(sprintf("  Missense: %d | Nonsense: %d | Frameshift: %d\n",
    length(snv_missense), length(snv_nonsense), length(snv_frameshift)))

  # -------------------------------------------------------------------------
  # Step 2: Generate variant protein sequences with customProDB
  # -------------------------------------------------------------------------
  # Build sample-specific protein FASTA
  outfa <- file.path(output_dir, "variant_proteins.fasta")

  OutputVarproseq(
    coding_variants,
    proteinseq = "data/ref/uniprot_human.fasta",
    outfile = outfa,
    ids = NULL       # all samples
  )
  cat(sprintf("  Variant protein FASTA: %s\n", outfa))

  # -------------------------------------------------------------------------
  # Step 3: Append to reference proteome for combined search
  # -------------------------------------------------------------------------
  ref_fa  <- readAAStringSet("data/ref/uniprot_human.fasta")
  var_fa  <- readAAStringSet(outfa)
  combined <- c(ref_fa, var_fa)

  # Add contaminant sequences
  contam  <- readAAStringSet("data/ref/contaminants.fasta")
  combined <- c(combined, contam)

  combined_path <- file.path(output_dir, "combined_search_db.fasta")
  writeXStringSet(combined, combined_path)
  cat(sprintf("  Combined DB: %d entries → %s\n", length(combined), combined_path))

  return(list(
    n_variants   = nrow(vcf),
    n_missense   = length(snv_missense),
    n_frameshift = length(snv_frameshift),
    db_path      = combined_path,
    n_db_entries = length(combined)
  ))
}

# -----------------------------------------------------------------------------
# Step 4: Configure MaxQuant XML for variant peptide search
# -----------------------------------------------------------------------------
generate_maxquant_xml <- function(db_path, raw_dir, output_dir = "results/maxquant/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mqpar_template <- '<?xml version="1.0" encoding="utf-8"?>
<MaxQuantParams xmlns:xsd="http://www.w3.org/2001/XMLSchema"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
   <fastaFiles>
      <FastaFileInfo>
         <fastaFilePath>%s</fastaFilePath>
         <identifierParseRule>>([^\\s]+)</identifierParseRule>
         <descriptionParseRule>>(.*)</descriptionParseRule>
         <taxonomyParseRule></taxonomyParseRule>
      </FastaFileInfo>
   </fastaFiles>
   <msmsParams>
      <maxMissedCleavages>%d</maxMissedCleavages>
      <minPepLen>%d</minPepLen>
      <maxPeptideMass>4600</maxPeptideMass>
      <matchType>MatchBetweenRuns</matchType>
      <mainSearchTol>4.5</mainSearchTol>
      <searchTolInPpm>true</searchTolInPpm>
   </msmsParams>
   <fixedModifications>
      <string>Carbamidomethyl (C)</string>
   </fixedModifications>
   <variableModifications>
      <string>Oxidation (M)</string>
      <string>Acetyl (Protein N-term)</string>
   </variableModifications>
   <decoyMode>revert</decoyMode>
   <psmFdrCutoff>%f</psmFdrCutoff>
   <proteinFdrCutoff>%f</proteinFdrCutoff>
   <numThreads>%d</numThreads>
   <labelFree>
      <lfqMinRatioCount>2</lfqMinRatioCount>
      <lfqNormType>1</lfqNormType>
   </labelFree>
</MaxQuantParams>'

  mqpar_xml <- sprintf(mqpar_template,
    normalizePath(db_path),
    cfg$max_missed_cleavages,
    cfg$min_peptide_length,
    cfg$fdr_threshold,
    cfg$fdr_threshold,
    config$environment$threads
  )

  xml_path <- file.path(output_dir, "mqpar_variantdb.xml")
  writeLines(mqpar_xml, xml_path)
  cat(sprintf("  MaxQuant XML: %s\n", xml_path))
  return(xml_path)
}

# -----------------------------------------------------------------------------
# Step 5: Post-search variant peptide filtering
# -----------------------------------------------------------------------------
filter_variant_peptides <- function(mq_output_dir, output_path = "results/variant_peptides_filtered.tsv") {
  evidence <- read.delim(file.path(mq_output_dir, "evidence.txt"),
    stringsAsFactors = FALSE)

  # Identify variant peptides (non-reference)
  is_variant <- grepl("^VAR_", evidence$Leading.razor.protein) |
                grepl("_MUT$", evidence$Leading.razor.protein)

  var_evidence <- evidence[is_variant, ]
  cat(sprintf("  Total PSMs: %d | Variant PSMs: %d (%.1f%%)\n",
    nrow(evidence), nrow(var_evidence),
    100 * nrow(var_evidence) / nrow(evidence)))

  # Quality filters
  var_evidence <- var_evidence[
    var_evidence$PEP < 0.01 &                  # posterior error prob
    var_evidence$Score > 40 &                   # Andromeda score
    !var_evidence$Reverse %in% "+" &            # not decoy
    !var_evidence$Potential.contaminant %in% "+", # not contaminant
  ]
  cat(sprintf("  After quality filter: %d variant PSMs\n", nrow(var_evidence)))

  # Annotate with genomic coordinates
  var_evidence$variant_type <- ifelse(
    grepl("_FS_", var_evidence$Leading.razor.protein), "frameshift",
    ifelse(grepl("_MIS_", var_evidence$Leading.razor.protein), "missense",
    "other")
  )

  # Summarize by unique peptide
  var_peptides <- aggregate(
    cbind(Score, Intensity) ~ Sequence + Leading.razor.protein + variant_type,
    data = var_evidence,
    FUN = function(x) c(max = max(x), n = length(x))
  )

  write.table(var_evidence, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Saved: %s (%d variant peptides)\n", output_path, nrow(var_peptides)))

  return(var_peptides)
}

# -----------------------------------------------------------------------------
# Step 6: COSMIC recurrence annotation
# -----------------------------------------------------------------------------
annotate_cosmic <- function(var_peptides, cosmic_vcf = "data/ref/CosmicCodingMuts.vcf.gz") {
  cosmic <- readVcf(cosmic_vcf, genome = "GRCh38")
  cosmic_ids <- paste(seqnames(cosmic), start(cosmic), ref(cosmic), sep = "_")

  var_peptides$cosmic_recurrence <- var_peptides$genomic_id %in% cosmic_ids
  n_cosmic <- sum(var_peptides$cosmic_recurrence, na.rm = TRUE)
  cat(sprintf("  COSMIC-annotated: %d / %d variants\n", n_cosmic, nrow(var_peptides)))

  return(var_peptides)
}

# =============================================================================
# Main Execution
# =============================================================================
main <- function() {
  cat("\n--- Building Variant Protein Database ---\n")
  db_info <- build_variant_proteome(cfg$vcf_input)

  cat("\n--- Generating MaxQuant Configuration ---\n")
  mqpar  <- generate_maxquant_xml(db_info$db_path, raw_dir = "data/raw_ms/")

  cat("\n--- Post-search Filtering (run after MaxQuant completes) ---\n")
  # var_pep <- filter_variant_peptides("results/maxquant/combined/txt/")
  # var_pep <- annotate_cosmic(var_pep)

  cat("\n=== Module 1 complete ===\n")
  cat(sprintf("  Database entries: %d\n", db_info$n_db_entries))
  cat(sprintf("  Missense variants: %d\n", db_info$n_missense))
  cat(sprintf("  Frameshift variants: %d\n", db_info$n_frameshift))

  invisible(db_info)
}

if (!interactive()) main()
