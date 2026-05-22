# ============================================================
# Step 1: Quality Control — adapter removal, dedup, host removal
# ============================================================

rule fastp_trim:
    """Adapter removal + quality trimming with fastp"""
    input:
        r1=lambda wc: samples_df.loc[samples_df["sample"]==wc.sample, "fq1"].values[0],
        r2=lambda wc: samples_df.loc[samples_df["sample"]==wc.sample, "fq2"].values[0],
    output:
        r1="results/qc/trimmed/{sample}_R1.fastq.gz",
        r2="results/qc/trimmed/{sample}_R2.fastq.gz",
        json="results/qc/fastp/{sample}_fastp.json",
        html="results/qc/fastp/{sample}_fastp.html",
    params:
        min_len=config["qc"]["min_length"],
        min_qual=config["qc"]["min_quality"],
        adapter=config["qc"]["adapter_file"],
    threads: config["resources"]["threads_qc"]
    conda: "../envs/qc.yaml"
    log: "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --json {output.json} --html {output.html} \
            --adapter_fasta {params.adapter} \
            --qualified_quality_phred {params.min_qual} \
            --length_required {params.min_len} \
            --low_complexity_filter \
            --complexity_threshold {config[qc][complexity_threshold]} \
            --thread {threads} \
            --detect_adapter_for_pe \
            2> {log}
        """

rule clumpify_dedup:
    """Optical + PCR duplicate removal with Clumpify (BBTools)"""
    input:
        r1="results/qc/trimmed/{sample}_R1.fastq.gz",
        r2="results/qc/trimmed/{sample}_R2.fastq.gz",
    output:
        r1="results/qc/dedup/{sample}_R1.fastq.gz",
        r2="results/qc/dedup/{sample}_R2.fastq.gz",
    params:
        dist=config["qc"]["dedup_optical_dist"],
    threads: config["resources"]["threads_qc"]
    conda: "../envs/qc.yaml"
    log: "logs/clumpify/{sample}.log"
    shell:
        """
        clumpify.sh \
            in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            dedupe=t optical=t dupedist={params.dist} \
            subs=0 passes=2 \
            2> {log}
        """

rule bowtie2_host_removal:
    """Remove host (human) reads via Bowtie2 mapping"""
    input:
        r1="results/qc/dedup/{sample}_R1.fastq.gz",
        r2="results/qc/dedup/{sample}_R2.fastq.gz",
    output:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
        stats="results/qc/host_stats/{sample}_host_mapping.txt",
    params:
        idx=config["databases"]["host_genome"],
        mapq=config["qc"]["host_map_quality"],
    threads: config["resources"]["threads_qc"]
    conda: "../envs/qc.yaml"
    log: "logs/bowtie2_host/{sample}.log"
    shell:
        """
        bowtie2 -x {params.idx} \
            -1 {input.r1} -2 {input.r2} \
            --very-sensitive -p {threads} \
            --un-conc-gz results/qc/hostfree/{wildcards.sample}_R%.fastq.gz \
            2> {output.stats} | \
        samtools view -bS -q {params.mapq} - > /dev/null 2>> {log}
        """

rule fastqc_clean:
    """FastQC on host-free reads"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
    output:
        directory("results/qc/fastqc/{sample}"),
    threads: 2
    conda: "../envs/qc.yaml"
    log: "logs/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {output}
        fastqc -t {threads} -o {output} {input.r1} {input.r2} 2> {log}
        """

rule multiqc:
    """Aggregate QC reports"""
    input:
        fastp=expand("results/qc/fastp/{sample}_fastp.json", sample=SAMPLES),
        fastqc=expand("results/qc/fastqc/{sample}", sample=SAMPLES),
        host=expand("results/qc/host_stats/{sample}_host_mapping.txt", sample=SAMPLES),
    output:
        "results/qc/multiqc_report.html",
    conda: "../envs/qc.yaml"
    log: "logs/multiqc.log"
    shell:
        """
        multiqc results/qc/ -o results/qc/ -f --no-data-dir 2> {log}
        """
