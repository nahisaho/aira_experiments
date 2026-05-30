# =============================================================================
# Taxonomic Classification Rules
# Kraken2, Bracken, MetaPhlAn4 comparison
# =============================================================================

rule kraken2_classify:
    """Taxonomic classification with Kraken2"""
    input:
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
    output:
        report="results/taxonomy/kraken2/{sample}.kreport",
        output="results/taxonomy/kraken2/{sample}.kraken2",
    params:
        db=config["databases"]["kraken2_db"],
        confidence=config["taxonomy"]["kraken2"]["confidence"],
        min_hit=config["taxonomy"]["kraken2"]["min_hit_groups"],
    threads: 8
    log:
        "logs/taxonomy/kraken2_{sample}.log"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --paired {input.r1} {input.r2} \
            --output {output.output} \
            --report {output.report} \
            --confidence {params.confidence} \
            --minimum-hit-groups {params.min_hit} \
            --threads {threads} \
            --gzip-compressed \
            2> {log}
        """

rule bracken_abundance:
    """Species-level abundance estimation with Bracken"""
    input:
        kreport="results/taxonomy/kraken2/{sample}.kreport",
    output:
        "results/taxonomy/bracken/{sample}.bracken",
    params:
        db=config["databases"]["kraken2_db"],
        threshold=config["taxonomy"]["bracken"]["threshold"],
        read_len=config["taxonomy"]["bracken"]["read_length"],
        level=config["taxonomy"]["bracken"]["level"],
    log:
        "logs/taxonomy/bracken_{sample}.log"
    shell:
        """
        bracken \
            -d {params.db} \
            -i {input.kreport} \
            -o {output} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} \
            2> {log}
        """

rule metaphlan4_profile:
    """Taxonomic profiling with MetaPhlAn4"""
    input:
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
    output:
        profile="results/taxonomy/metaphlan4/{sample}_profile.tsv",
        bowtie2out="results/taxonomy/metaphlan4/{sample}.bowtie2.bz2",
    params:
        db=config["databases"]["metaphlan_db"],
        analysis=config["taxonomy"]["metaphlan4"]["analysis_type"],
        min_mapq=config["taxonomy"]["metaphlan4"]["min_mapq"],
    threads: 8
    log:
        "logs/taxonomy/metaphlan4_{sample}.log"
    shell:
        """
        metaphlan \
            {input.r1},{input.r2} \
            --input_type fastq \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2out} \
            --output_file {output.profile} \
            -t {params.analysis} \
            --min_mapq_val {params.min_mapq} \
            --nproc {threads} \
            2> {log}
        """

rule merge_metaphlan:
    """Merge MetaPhlAn4 profiles"""
    input:
        expand("results/taxonomy/metaphlan4/{sample}_profile.tsv", sample=SAMPLE_IDS),
    output:
        "results/taxonomy/merged_metaphlan_profiles.tsv",
    log:
        "logs/taxonomy/merge_metaphlan.log"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output} 2> {log}
        """

rule merge_kraken2:
    """Merge Kraken2/Bracken profiles into a unified table"""
    input:
        expand("results/taxonomy/bracken/{sample}.bracken", sample=SAMPLE_IDS),
    output:
        "results/taxonomy/merged_kraken2_profiles.tsv",
    log:
        "logs/taxonomy/merge_kraken2.log"
    script:
        "../scripts/merge_kraken2.py"
