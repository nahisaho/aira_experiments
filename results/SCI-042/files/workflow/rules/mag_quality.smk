# ============================================================
# Step 5: MAG Quality Assessment & Phylogenetic Placement
#         CheckM2 + GTDB-Tk + QUAST
# ============================================================

rule checkm2_evaluate:
    """CheckM2: ML-based completeness & contamination estimation"""
    input:
        bins="results/binning/dastool/all_bins",
    output:
        "results/mag_quality/checkm2_report.tsv",
    params:
        db=config["databases"]["checkm2_db"],
        outdir="results/mag_quality/checkm2",
    threads: config["resources"]["threads_annotation"]
    conda: "../envs/mag_quality.yaml"
    log: "logs/checkm2.log"
    shell:
        """
        checkm2 predict \
            --input {input.bins} \
            --output-directory {params.outdir} \
            --database_path {params.db} \
            --threads {threads} \
            --extension fa \
            --force \
            2> {log}
        cp {params.outdir}/quality_report.tsv {output}
        """

rule filter_quality_mags:
    """Filter MAGs by MIMAG medium-quality thresholds"""
    input:
        report="results/mag_quality/checkm2_report.tsv",
        bins="results/binning/dastool/all_bins",
    output:
        directory("results/mag_quality/filtered_mags"),
        "results/mag_quality/mag_summary.tsv",
    params:
        comp=config["mag_quality"]["completeness_threshold"],
        cont=config["mag_quality"]["contamination_threshold"],
        hq_comp=config["mag_quality"]["high_quality_comp"],
        hq_cont=config["mag_quality"]["high_quality_cont"],
    conda: "../envs/mag_quality.yaml"
    log: "logs/filter_mags.log"
    script:
        "../scripts/filter_mags.py"

rule gtdbtk_classify:
    """GTDB-Tk: phylogenetic placement of filtered MAGs"""
    input:
        mags="results/mag_quality/filtered_mags",
    output:
        directory("results/mag_quality/gtdbtk_output"),
    params:
        db=config["databases"]["gtdbtk_db"],
    threads: config["resources"]["threads_annotation"]
    resources:
        mem_mb=config["resources"]["mem_mb_binning"],
    conda: "../envs/mag_quality.yaml"
    log: "logs/gtdbtk.log"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf \
            --genome_dir {input.mags} \
            --out_dir {output} \
            --extension fa \
            --cpus {threads} \
            --pplacer_cpus 1 \
            --skip_ani_screen \
            2> {log}
        """
