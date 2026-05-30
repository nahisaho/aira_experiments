# =============================================================================
# Assembly Rules
# MEGAHIT metagenomic assembly
# =============================================================================

rule megahit_assembly:
    """Metagenomic assembly with MEGAHIT"""
    input:
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
    output:
        contigs="results/assembly/{sample}/final.contigs.fa",
    params:
        outdir="results/assembly/{sample}",
        min_contig=config["assembly"]["min_contig_length"],
        k_list=config["assembly"]["k_list"],
        preset=config["assembly"]["megahit_preset"],
    threads: 16
    log:
        "logs/assembly/megahit_{sample}.log"
    shell:
        """
        rm -rf {params.outdir}
        megahit \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            --min-contig-len {params.min_contig} \
            --k-list {params.k_list} \
            --presets {params.preset} \
            --num-cpu-threads {threads} \
            2> {log}
        """

rule contig_stats:
    """Calculate assembly statistics"""
    input:
        "results/assembly/{sample}/final.contigs.fa",
    output:
        "results/assembly/{sample}/assembly_stats.tsv",
    log:
        "logs/assembly/stats_{sample}.log"
    script:
        "../scripts/assembly_stats.py"
