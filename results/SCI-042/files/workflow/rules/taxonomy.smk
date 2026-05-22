# ============================================================
# Step 2: Assembly-Free Taxonomic Classification
#         Kraken2 + Bracken vs MetaPhlAn 4 — comparative profiling
# ============================================================

rule kraken2_classify:
    """Kraken2 k-mer based classification"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
    output:
        report="results/taxonomy/kraken2/{sample}.kreport",
        output="results/taxonomy/kraken2/{sample}.kraken2",
    params:
        db=config["databases"]["kraken2_db"],
        confidence=0.2,
    threads: config["resources"]["threads_qc"]
    conda: "../envs/taxonomy.yaml"
    log: "logs/kraken2/{sample}.log"
    shell:
        """
        kraken2 --db {params.db} \
            --paired {input.r1} {input.r2} \
            --output {output.output} \
            --report {output.report} \
            --confidence {params.confidence} \
            --threads {threads} \
            --gzip-compressed \
            --memory-mapping \
            2> {log}
        """

rule bracken_reestimation:
    """Bracken: Bayesian re-estimation of Kraken2 abundances"""
    input:
        kreport="results/taxonomy/kraken2/{sample}.kreport",
    output:
        bracken="results/taxonomy/bracken/{sample}.bracken",
        kreport="results/taxonomy/bracken/{sample}_bracken.kreport",
    params:
        db=config["databases"]["kraken2_db"],
        level="S",       # Species level
        threshold=10,     # Minimum reads
    conda: "../envs/taxonomy.yaml"
    log: "logs/bracken/{sample}.log"
    shell:
        """
        bracken -d {params.db} \
            -i {input.kreport} \
            -o {output.bracken} \
            -w {output.kreport} \
            -r 150 -l {params.level} -t {params.threshold} \
            2> {log}
        """

rule metaphlan4_profile:
    """MetaPhlAn 4: clade-specific marker gene profiling"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
    output:
        profile="results/taxonomy/metaphlan4/{sample}_profile.tsv",
        bowtie2="results/taxonomy/metaphlan4/{sample}_bowtie2.bz2",
    params:
        db=config["databases"]["metaphlan_db"],
    threads: config["resources"]["threads_qc"]
    conda: "../envs/taxonomy.yaml"
    log: "logs/metaphlan4/{sample}.log"
    shell:
        """
        metaphlan {input.r1},{input.r2} \
            --input_type fastq \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2} \
            --nproc {threads} \
            --unclassified_estimation \
            --tax_lev a \
            -o {output.profile} \
            2> {log}
        """

rule merge_metaphlan4:
    """Merge MetaPhlAn 4 profiles across samples"""
    input:
        expand("results/taxonomy/metaphlan4/{sample}_profile.tsv", sample=SAMPLES),
    output:
        "results/taxonomy/merged_metaphlan4.tsv",
    conda: "../envs/taxonomy.yaml"
    log: "logs/merge_metaphlan4.log"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output} 2> {log}
        """

rule taxonomy_comparison:
    """Compare Kraken2/Bracken vs MetaPhlAn4 concordance"""
    input:
        bracken=expand("results/taxonomy/bracken/{sample}.bracken", sample=SAMPLES),
        metaphlan=expand("results/taxonomy/metaphlan4/{sample}_profile.tsv", sample=SAMPLES),
    output:
        "results/taxonomy/classifier_comparison.tsv",
    conda: "../envs/taxonomy.yaml"
    log: "logs/taxonomy_comparison.log"
    script:
        "../scripts/compare_classifiers.py"
