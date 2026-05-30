# =============================================================================
# Quality Control Rules
# Host removal, adapter trimming, deduplication
# =============================================================================

rule fastp_trim:
    """Adapter trimming and quality filtering with fastp"""
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["r1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["r2"],
    output:
        r1="results/qc/{sample}_trimmed_R1.fastq.gz",
        r2="results/qc/{sample}_trimmed_R2.fastq.gz",
        json="results/qc/{sample}_fastp.json",
        html="results/qc/{sample}_fastp.html",
    params:
        quality=config["qc"]["min_quality"],
        length=config["qc"]["min_length"],
        front=config["qc"]["trim_front"],
        tail=config["qc"]["trim_tail"],
        adapter=config["qc"]["adapter_file"],
    threads: 4
    log:
        "logs/qc/fastp_{sample}.log"
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --json {output.json} --html {output.html} \
            --qualified_quality_phred {params.quality} \
            --length_required {params.length} \
            --trim_front1 {params.front} --trim_front2 {params.front} \
            --trim_tail1 {params.tail} --trim_tail2 {params.tail} \
            --adapter_fasta {params.adapter} \
            --dedup --dup_calc_accuracy {config[qc][dedup_accuracy]} \
            --thread {threads} \
            2> {log}
        """

rule host_removal_index:
    """Build Bowtie2 index for host genome"""
    input:
        config["databases"]["host_genome"] + ".fa"
    output:
        touch("results/qc/.host_index_done")
    params:
        prefix=config["databases"]["host_genome"]
    threads: 8
    log:
        "logs/qc/bowtie2_build.log"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input} {params.prefix} \
            2> {log}
        """

rule host_removal:
    """Remove host (human) reads using Bowtie2"""
    input:
        r1="results/qc/{sample}_trimmed_R1.fastq.gz",
        r2="results/qc/{sample}_trimmed_R2.fastq.gz",
        index="results/qc/.host_index_done",
    output:
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
        host_bam="results/qc/{sample}_host_aligned.bam",
    params:
        index=config["databases"]["host_genome"],
        sensitivity=config["qc"]["host_removal_sensitivity"],
    threads: 8
    log:
        "logs/qc/host_removal_{sample}.log"
    shell:
        """
        bowtie2 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            {params.sensitivity} \
            --threads {threads} \
            --un-conc-gz results/qc/{wildcards.sample}_R%.clean.fastq.gz \
            2> {log} \
        | samtools view -bS -f 4 - > {output.host_bam}
        """

rule multiqc:
    """Aggregate QC reports"""
    input:
        expand("results/qc/{sample}_fastp.json", sample=SAMPLE_IDS),
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/qc/multiqc.log"
    shell:
        """
        multiqc results/qc/ -o results/qc/ --force 2> {log}
        """
