# =============================================================================
# MAG Quality Assessment and Phylogenetic Placement
# CheckM2 and GTDB-Tk
# =============================================================================

rule checkm2_assess:
    """Assess MAG quality with CheckM2"""
    input:
        bins="results/binning/dastool/{sample}_DASTool_bins",
    output:
        "results/mags/{sample}/checkm2_results.tsv",
    params:
        db=config["databases"]["checkm2_db"],
        outdir="results/mags/{sample}/checkm2",
    threads: config["mag_quality"]["checkm2"]["threads"]
    log:
        "logs/mags/checkm2_{sample}.log"
    shell:
        """
        checkm2 predict \
            --input {input.bins} \
            --output-directory {params.outdir} \
            --database_path {params.db} \
            --threads {threads} \
            -x fa \
            2> {log}
        cp {params.outdir}/quality_report.tsv {output}
        """

rule gtdbtk_classify:
    """Taxonomic classification of MAGs with GTDB-Tk"""
    input:
        bins="results/binning/dastool/{sample}_DASTool_bins",
    output:
        "results/mags/{sample}/gtdbtk_results.tsv",
    params:
        db=config["databases"]["gtdbtk_db"],
        outdir="results/mags/{sample}/gtdbtk",
    threads: config["mag_quality"]["gtdbtk"]["threads"]
    log:
        "logs/mags/gtdbtk_{sample}.log"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf \
            --genome_dir {input.bins} \
            --out_dir {params.outdir} \
            --extension fa \
            --cpus {threads} \
            2> {log}
        # Merge bacterial and archaeal classifications
        cat {params.outdir}/classify/gtdbtk.*.summary.tsv | head -1 > {output}
        cat {params.outdir}/classify/gtdbtk.*.summary.tsv | grep -v "^user_genome" >> {output} || true
        """

rule filter_mags:
    """Filter MAGs by quality thresholds and summarize"""
    input:
        checkm=expand("results/mags/{sample}/checkm2_results.tsv", sample=SAMPLE_IDS),
        gtdbtk=expand("results/mags/{sample}/gtdbtk_results.tsv", sample=SAMPLE_IDS),
    output:
        "results/mags/all_mags_summary.tsv",
    params:
        min_completeness=config["mag_quality"]["completeness_threshold"],
        max_contamination=config["mag_quality"]["contamination_threshold"],
        hq_completeness=config["mag_quality"]["high_quality_completeness"],
        hq_contamination=config["mag_quality"]["high_quality_contamination"],
    log:
        "logs/mags/filter_mags.log"
    script:
        "../scripts/filter_mags.py"
