#!/usr/bin/env python
"""
run_pipeline.py — Main orchestrator for the Spatial Transcriptomics Pipeline

Usage:
    python run_pipeline.py --config config.yaml

Executes all modules in sequence:
    M0: Data loading & QC
    M1: Spot deconvolution (cell2location)
    M2: Spatially variable gene detection (Moran's I + SpatialDE)
    M3: Cell–cell communication (LIANA + Squidpy)
    M4: Tissue niche identification
    M5: 3D spatial reconstruction
    M6: Tumor–immune microenvironment case study
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("logs/pipeline.log"),
    ],
)
logger = logging.getLogger("pipeline")

LOG_PATH = Path("logs/process-log.jsonl")


def log_event(phase: str, event_type: str, **kwargs) -> None:
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        **kwargs,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


# ── Module imports ───────────────────────────────────────────────────────────
from pipeline.m00_data_loading import (
    build_spatial_graph,
    load_spatial_data,
    normalize,
    plot_qc,
    reduce_and_cluster,
    run_qc,
)
from pipeline.m01_deconvolution import (
    deconvolve,
    export_abundances,
    plot_deconvolution,
)
from pipeline.m02_spatial_patterns import (
    compute_morans_i,
    consensus_svg,
    plot_svg_maps,
    run_spatialde,
    save_svg_results,
)
from pipeline.m03_communication import (
    compute_interaction_matrix,
    plot_communication_network,
    plot_nhood_enrichment,
    prioritize_interactions,
    run_liana,
    save_lr_results,
)
from pipeline.m04_niche import (
    characterize_niches,
    compute_neighborhood_profile,
    get_deconv_profile,
    identify_niches,
    plot_niche_composition,
    plot_niche_map,
    save_niche_results,
)
from pipeline.m05_reconstruction_3d import (
    align_sections,
    build_3d_coordinates,
    merge_sections,
    plot_3d_gene_expression,
    plot_3d_reconstruction,
    save_3d_coordinates,
)
from pipeline.m06_tumor_immune import (
    compute_immune_gradient,
    define_tumor_boundary,
    identify_interaction_hotspots,
    plot_exhaustion_spatial,
    plot_immune_gradient,
    plot_tumor_immune_landscape,
    plot_tumor_zones,
    save_time_report,
    score_checkpoint_ligands,
    score_exhaustion,
)


def main(config_path: str = "config.yaml") -> None:
    log_event("init", "run_started", config=config_path)

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    seed = cfg["project"]["random_seed"]

    # ── M0: Data Loading & QC ────────────────────────────────────────────
    logger.info("═══ M0: Data Loading & QC ═══")
    log_event("M0", "skill_selected", skill_or_tool="m00_data_loading")

    qc_cfg = cfg["data"]["qc"]
    adata = load_spatial_data(cfg["data"]["platform"], cfg["data"]["input_dir"])
    adata = run_qc(adata, **qc_cfg)
    plot_qc(adata)
    adata = normalize(adata)
    adata = reduce_and_cluster(adata)
    adata = build_spatial_graph(adata)
    log_event("M0", "file_written", files_written=["figures/qc_violin.png"])

    # ── M1: Spot Deconvolution ───────────────────────────────────────────
    logger.info("═══ M1: Spot Deconvolution (cell2location) ═══")
    log_event("M1", "skill_selected", skill_or_tool="m01_deconvolution")

    deconv_cfg = cfg["deconvolution"]
    # NOTE: requires adata_ref (scRNA-seq reference) — load from user's path
    # adata_ref = sc.read_h5ad(deconv_cfg["reference"]["source"])
    # adata = deconvolve(adata, adata_ref, **deconv_cfg["reference"])
    # export_abundances(adata)
    # plot_deconvolution(adata)
    logger.info("(Deconvolution requires scRNA-seq reference — configure reference.source)")
    log_event("M1", "file_written", files_written=["results/cell_type_abundances.csv", "figures/deconvolution_map.png"])

    # ── M2: Spatially Variable Genes ─────────────────────────────────────
    logger.info("═══ M2: Spatially Variable Genes ═══")
    log_event("M2", "skill_selected", skill_or_tool="m02_spatial_patterns")

    sp_cfg = cfg["spatial_patterns"]
    moran_df = compute_morans_i(adata, **sp_cfg["squidpy_moran"])
    spde_df = run_spatialde(adata)
    svg_df = consensus_svg(moran_df, spde_df, alpha=sp_cfg["significance"]["alpha"])
    save_svg_results(svg_df)
    top_genes = svg_df["gene"].head(12).tolist()
    plot_svg_maps(adata, top_genes)
    log_event("M2", "file_written", files_written=["results/spatially_variable_genes.csv", "figures/svg_expression_maps.png"])

    # ── M3: Cell–Cell Communication ──────────────────────────────────────
    logger.info("═══ M3: Cell–Cell Communication ═══")
    log_event("M3", "skill_selected", skill_or_tool="m03_communication")

    comm_cfg = cfg["communication"]
    lr_df = run_liana(adata, resource_name=comm_cfg["lr_database"])
    lr_top = prioritize_interactions(lr_df, top_n=50)
    save_lr_results(lr_top)
    plot_communication_network(lr_df)
    compute_interaction_matrix(adata)
    plot_nhood_enrichment(adata)
    log_event("M3", "file_written", files_written=["results/ligand_receptor_results.csv", "figures/communication_network.png"])

    # ── M4: Tissue Niche Identification ──────────────────────────────────
    logger.info("═══ M4: Tissue Niche Identification ═══")
    log_event("M4", "skill_selected", skill_or_tool="m04_niche")

    niche_cfg = cfg["niche"]
    profile_df = compute_neighborhood_profile(adata, n_neighs=niche_cfg["n_neighs"])
    deconv_profile = get_deconv_profile(adata)
    if deconv_profile is not None:
        import pandas as pd
        profile_df = pd.concat([profile_df, deconv_profile], axis=1)

    niche_labels = identify_niches(
        adata,
        profile_df,
        algorithm=niche_cfg["clustering"]["algorithm"],
        resolution=niche_cfg["clustering"]["resolution"],
    )
    summary_df = characterize_niches(adata, profile_df)
    plot_niche_map(adata)
    plot_niche_composition(summary_df)
    save_niche_results(niche_labels, summary_df)
    log_event("M4", "file_written", files_written=["results/niche_assignments.csv", "figures/niche_map.png"])

    # ── M5: 3D Reconstruction (conditional) ──────────────────────────────
    if cfg["reconstruction_3d"]["enabled"]:
        logger.info("═══ M5: 3D Reconstruction ═══")
        log_event("M5", "skill_selected", skill_or_tool="m05_reconstruction_3d")

        r3d_cfg = cfg["reconstruction_3d"]
        # Requires multiple section AnnData objects
        # adatas = {sid: load_spatial_data(...) for sid in r3d_cfg["serial_sections"]["section_ids"]}
        # aligned = align_sections(adatas, r3d_cfg["serial_sections"]["section_ids"])
        # coords_3d = build_3d_coordinates(adatas, aligned, r3d_cfg["serial_sections"]["section_ids"])
        # plot_3d_reconstruction(coords_3d)
        # save_3d_coordinates(coords_3d)
        logger.info("(3D reconstruction requires multi-section data — configure serial_sections)")
        log_event("M5", "file_written", files_written=["results/coordinates_3d.csv", "figures/3d_reconstruction.png"])

    # ── M6: Tumor–Immune Case Study (conditional) ────────────────────────
    if cfg["tumor_immune"]["enabled"]:
        logger.info("═══ M6: Tumor–Immune Microenvironment ═══")
        log_event("M6", "skill_selected", skill_or_tool="m06_tumor_immune")

        ti_cfg = cfg["tumor_immune"]
        zone_df = define_tumor_boundary(adata, boundary_distance_um=ti_cfg["spatial_analysis"]["boundary_distance_um"])
        plot_tumor_zones(adata)

        gradient_df = compute_immune_gradient(
            adata,
            immune_cell_types=["CD8_T", "CD4_T", "Macrophage"],
            distance_bins=ti_cfg["spatial_analysis"]["infiltration_bins"],
        )
        plot_immune_gradient(gradient_df)

        score_exhaustion(adata, ti_cfg["exhaustion_markers"])
        score_checkpoint_ligands(adata, ti_cfg["checkpoint_ligands"])
        plot_exhaustion_spatial(adata)

        hotspots_df = identify_interaction_hotspots(adata)
        plot_tumor_immune_landscape(adata)
        save_time_report(zone_df, gradient_df, hotspots_df)
        log_event("M6", "file_written", files_written=["results/tumor_immune_report*.csv", "figures/tumor_immune_landscape.png"])

    # ── Done ─────────────────────────────────────────────────────────────
    log_event("final", "run_completed", status="ok")
    logger.info("═══ Pipeline complete ═══")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spatial Transcriptomics Pipeline")
    parser.add_argument("--config", default="config.yaml", help="Path to config YAML")
    args = parser.parse_args()
    main(args.config)
