#!/usr/bin/env bash
# =============================================================================
# Cancer Proteogenomics Integrated Analysis Pipeline — Master Runner
# =============================================================================
# Usage: bash scripts/00_run_pipeline.sh [module_number]
#   e.g. bash scripts/00_run_pipeline.sh       # run all
#        bash scripts/00_run_pipeline.sh 3      # run only module 3
# =============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"
LOG_DIR="${PROJECT_DIR}/logs"
RESULTS_DIR="${PROJECT_DIR}/results"
FIGURES_DIR="${PROJECT_DIR}/figures"

mkdir -p "$LOG_DIR" "$RESULTS_DIR" "$FIGURES_DIR"

echo "============================================================"
echo " Cancer Proteogenomics Pipeline"
echo " Start: $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project: ${PROJECT_DIR}"
echo "============================================================"

run_module() {
    local num=$1
    local script=$2
    local desc=$3
    local ext="${script##*.}"

    echo ""
    echo "--- Module ${num}: ${desc} ---"
    echo "  Script: ${script}"
    echo "  Start:  $(date '+%H:%M:%S')"

    local logfile="${LOG_DIR}/module${num}_$(date '+%Y%m%d_%H%M%S').log"

    if [ "$ext" = "R" ]; then
        Rscript "$script" 2>&1 | tee "$logfile"
    elif [ "$ext" = "py" ]; then
        python3 "$script" 2>&1 | tee "$logfile"
    else
        bash "$script" 2>&1 | tee "$logfile"
    fi

    echo "  End:    $(date '+%H:%M:%S')"
    echo "  Log:    ${logfile}"
}

MODULE="${1:-all}"

if [ "$MODULE" = "all" ] || [ "$MODULE" = "1" ]; then
    run_module 1 "${PIPELINE_DIR}/01_variant_peptide_search.R" \
        "Variant Peptide Database & Search"
fi

if [ "$MODULE" = "all" ] || [ "$MODULE" = "2" ]; then
    run_module 2 "${PIPELINE_DIR}/02_rna_protein_discordance.R" \
        "RNA–Protein Discordance Analysis"
fi

if [ "$MODULE" = "all" ] || [ "$MODULE" = "3" ]; then
    run_module 3 "${PIPELINE_DIR}/03_phosphoproteomics_kinase.R" \
        "Phosphoproteomics & Kinase Activity"
fi

if [ "$MODULE" = "all" ] || [ "$MODULE" = "4" ]; then
    run_module 4 "${PIPELINE_DIR}/04_neoantigen_verification.py" \
        "Neoantigen Proteomics Verification"
fi

if [ "$MODULE" = "all" ] || [ "$MODULE" = "5" ]; then
    run_module 5 "${PIPELINE_DIR}/05_mofa_integration.R" \
        "MOFA+ Multi-Omics Integration"
fi

if [ "$MODULE" = "all" ] || [ "$MODULE" = "6" ]; then
    run_module 6 "${PIPELINE_DIR}/06_cptac_pdac_casestudy.R" \
        "CPTAC PDAC Case Study"
fi

echo ""
echo "============================================================"
echo " Pipeline Complete: $(date '+%Y-%m-%d %H:%M:%S')"
echo " Results: ${RESULTS_DIR}/"
echo " Figures: ${FIGURES_DIR}/"
echo " Logs:    ${LOG_DIR}/"
echo "============================================================"
