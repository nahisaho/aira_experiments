# =============================================================================
# Functional Annotation Rules
# HUMAnN3 and eggNOG-mapper integration
# =============================================================================

rule humann3_profile:
    """Functional profiling with HUMAnN3"""
    input:
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
        profile="results/taxonomy/metaphlan4/{sample}_profile.tsv",
    output:
        genefamilies="results/functional/humann3/{sample}_genefamilies.tsv",
        pathabundance="results/functional/humann3/{sample}_pathabundance.tsv",
        pathcoverage="results/functional/humann3/{sample}_pathcoverage.tsv",
    params:
        outdir="results/functional/humann3/",
        search=config["functional"]["humann3"]["search_mode"],
        db=config["databases"]["humann3_db"],
    threads: config["functional"]["humann3"]["threads"]
    log:
        "logs/functional/humann3_{sample}.log"
    shell:
        """
        # Concatenate paired reads for HUMAnN3
        cat {input.r1} {input.r2} > results/functional/humann3/{wildcards.sample}_concat.fastq.gz

        humann \
            --input results/functional/humann3/{wildcards.sample}_concat.fastq.gz \
            --output {params.outdir} \
            --taxonomic-profile {input.profile} \
            --search-mode {params.search} \
            --nucleotide-database {params.db}/chocophlan \
            --protein-database {params.db}/uniref \
            --threads {threads} \
            2> {log}

        rm -f results/functional/humann3/{wildcards.sample}_concat.fastq.gz
        """

rule humann3_normalize:
    """Normalize HUMAnN3 gene family and pathway abundance tables"""
    input:
        genefamilies="results/functional/humann3/{sample}_genefamilies.tsv",
        pathabundance="results/functional/humann3/{sample}_pathabundance.tsv",
    output:
        genefamilies_cpm="results/functional/humann3/{sample}_genefamilies_cpm.tsv",
        pathabundance_relab="results/functional/humann3/{sample}_pathabundance_relab.tsv",
    log:
        "logs/functional/humann3_normalize_{sample}.log"
    shell:
        """
        humann_renorm_table \
            --input {input.genefamilies} \
            --output {output.genefamilies_cpm} \
            --units cpm \
            2> {log}

        humann_renorm_table \
            --input {input.pathabundance} \
            --output {output.pathabundance_relab} \
            --units relab \
            2>> {log}
        """

rule merge_humann3:
    """Merge HUMAnN3 output tables across samples"""
    input:
        gf=expand("results/functional/humann3/{sample}_genefamilies.tsv", sample=SAMPLE_IDS),
        pa=expand("results/functional/humann3/{sample}_pathabundance.tsv", sample=SAMPLE_IDS),
    output:
        gf="results/functional/merged_genefamilies.tsv",
        pa="results/functional/merged_pathabundance.tsv",
    log:
        "logs/functional/merge_humann3.log"
    shell:
        """
        humann_join_tables \
            --input results/functional/humann3/ \
            --output {output.gf} \
            --file_name genefamilies \
            2> {log}

        humann_join_tables \
            --input results/functional/humann3/ \
            --output {output.pa} \
            --file_name pathabundance \
            2>> {log}
        """

rule prodigal_predict:
    """Gene prediction from assembled contigs"""
    input:
        "results/assembly/{sample}/final.contigs.fa",
    output:
        proteins="results/functional/prodigal/{sample}_proteins.faa",
        genes="results/functional/prodigal/{sample}_genes.gff",
    log:
        "logs/functional/prodigal_{sample}.log"
    shell:
        """
        prodigal \
            -i {input} \
            -a {output.proteins} \
            -o {output.genes} \
            -p meta \
            -f gff \
            2> {log}
        """

rule eggnog_annotate:
    """Functional annotation with eggNOG-mapper"""
    input:
        "results/functional/prodigal/{sample}_proteins.faa",
    output:
        "results/functional/eggnog/{sample}.emapper.annotations",
    params:
        db=config["databases"]["eggnog_db"],
        prefix="results/functional/eggnog/{sample}",
        tax_scope=config["functional"]["eggnog"]["tax_scope"],
    threads: config["functional"]["eggnog"]["threads"]
    log:
        "logs/functional/eggnog_{sample}.log"
    shell:
        """
        emapper.py \
            -i {input} \
            --output {params.prefix} \
            --data_dir {params.db} \
            --cpu {threads} \
            --tax_scope {params.tax_scope} \
            -m diamond \
            --override \
            2> {log}
        """
