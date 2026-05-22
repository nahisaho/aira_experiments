# ============================================================
# Step 4: Co-assembly + Genome Binning
#         MEGAHIT → MetaBAT2 / CONCOCT / MaxBin2 → DAS Tool
# ============================================================

rule megahit_assembly:
    """MEGAHIT co-assembly per sample"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
    output:
        contigs="results/assembly/{sample}/final_contigs.fa",
    params:
        outdir="results/assembly/{sample}",
        min_len=config["assembly"]["megahit"]["min_contig_len"],
        k_list=config["assembly"]["megahit"]["k_list"],
    threads: config["resources"]["threads_assembly"]
    resources:
        mem_mb=config["resources"]["mem_mb_assembly"],
    conda: "../envs/assembly.yaml"
    log: "logs/megahit/{sample}.log"
    shell:
        """
        rm -rf {params.outdir}/megahit_tmp
        megahit -1 {input.r1} -2 {input.r2} \
            -o {params.outdir}/megahit_tmp \
            --k-list {params.k_list} \
            --min-contig-len {params.min_len} \
            -t {threads} \
            -m 0.9 \
            2> {log}
        mv {params.outdir}/megahit_tmp/final.contigs.fa {output.contigs}
        rm -rf {params.outdir}/megahit_tmp
        """

rule assembly_stats:
    """QUAST assembly quality assessment"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
    output:
        "results/assembly/{sample}/quast/report.tsv",
    threads: 4
    conda: "../envs/assembly.yaml"
    log: "logs/quast/{sample}.log"
    shell:
        """
        quast {input.contigs} \
            -o results/assembly/{wildcards.sample}/quast \
            --min-contig 1000 \
            -t {threads} \
            2> {log}
        """

rule bowtie2_map_contigs:
    """Map reads back to assembled contigs for coverage"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
        contigs="results/assembly/{sample}/final_contigs.fa",
    output:
        bam="results/binning/mapping/{sample}.sorted.bam",
        bai="results/binning/mapping/{sample}.sorted.bam.bai",
    threads: config["resources"]["threads_binning"]
    conda: "../envs/assembly.yaml"
    log: "logs/bowtie2_contigs/{sample}.log"
    shell:
        """
        bowtie2-build {input.contigs} {input.contigs} 2>> {log}
        bowtie2 -x {input.contigs} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --very-sensitive-local \
            2>> {log} | \
        samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index {output.bam} 2>> {log}
        """

rule metabat2_depth:
    """Calculate contig depth for MetaBAT2"""
    input:
        bam="results/binning/mapping/{sample}.sorted.bam",
    output:
        "results/binning/metabat2/{sample}_depth.txt",
    conda: "../envs/binning.yaml"
    log: "logs/metabat2_depth/{sample}.log"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output} {input.bam} 2> {log}
        """

rule metabat2_bin:
    """MetaBAT2 binning"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
        depth="results/binning/metabat2/{sample}_depth.txt",
    output:
        directory("results/binning/metabat2/{sample}_bins"),
    params:
        min_size=config["binning"]["metabat2"]["min_bin_size"],
        min_len=config["binning"]["min_contig_length"],
    threads: config["resources"]["threads_binning"]
    conda: "../envs/binning.yaml"
    log: "logs/metabat2/{sample}.log"
    shell:
        """
        mkdir -p {output}
        metabat2 -i {input.contigs} \
            -a {input.depth} \
            -o {output}/bin \
            -m {params.min_len} \
            --minClsSize {params.min_size} \
            -t {threads} \
            2> {log}
        """

rule concoct_bin:
    """CONCOCT binning"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
        bam="results/binning/mapping/{sample}.sorted.bam",
    output:
        directory("results/binning/concoct/{sample}_bins"),
    params:
        chunk=config["binning"]["concoct"]["chunk_size"],
    threads: config["resources"]["threads_binning"]
    conda: "../envs/binning.yaml"
    log: "logs/concoct/{sample}.log"
    shell:
        """
        mkdir -p results/binning/concoct/{wildcards.sample}_tmp {output}

        cut_up_fasta.py {input.contigs} \
            -c {params.chunk} -o 0 --merge_last \
            -b results/binning/concoct/{wildcards.sample}_tmp/contigs_10K.bed \
            > results/binning/concoct/{wildcards.sample}_tmp/contigs_10K.fa

        concoct_coverage_table.py \
            results/binning/concoct/{wildcards.sample}_tmp/contigs_10K.bed \
            {input.bam} \
            > results/binning/concoct/{wildcards.sample}_tmp/coverage.tsv

        concoct \
            --composition_file results/binning/concoct/{wildcards.sample}_tmp/contigs_10K.fa \
            --coverage_file results/binning/concoct/{wildcards.sample}_tmp/coverage.tsv \
            -t {threads} \
            -b results/binning/concoct/{wildcards.sample}_tmp/ \
            2> {log}

        merge_cutup_clustering.py \
            results/binning/concoct/{wildcards.sample}_tmp/clustering_gt1000.csv \
            > results/binning/concoct/{wildcards.sample}_tmp/clustering_merged.csv

        extract_fasta_bins.py {input.contigs} \
            results/binning/concoct/{wildcards.sample}_tmp/clustering_merged.csv \
            --output_path {output}

        rm -rf results/binning/concoct/{wildcards.sample}_tmp
        """

rule maxbin2_bin:
    """MaxBin2 binning"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
        depth="results/binning/metabat2/{sample}_depth.txt",
    output:
        directory("results/binning/maxbin2/{sample}_bins"),
    params:
        min_len=config["binning"]["maxbin2"]["min_contig_length"],
        prob=config["binning"]["maxbin2"]["prob_threshold"],
    threads: config["resources"]["threads_binning"]
    conda: "../envs/binning.yaml"
    log: "logs/maxbin2/{sample}.log"
    shell:
        """
        mkdir -p {output}
        # Extract abundance from MetaBAT2 depth file
        cut -f1,3 {input.depth} | tail -n+2 \
            > results/binning/maxbin2/{wildcards.sample}_abundance.txt

        run_MaxBin.pl \
            -contig {input.contigs} \
            -abund results/binning/maxbin2/{wildcards.sample}_abundance.txt \
            -out {output}/bin \
            -min_contig_length {params.min_len} \
            -prob_threshold {params.prob} \
            -thread {threads} \
            2> {log}
        """

rule dastool_refine:
    """DAS Tool: consensus binning from MetaBAT2 + CONCOCT + MaxBin2"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
        metabat2="results/binning/metabat2/{sample}_bins",
        concoct="results/binning/concoct/{sample}_bins",
        maxbin2="results/binning/maxbin2/{sample}_bins",
    output:
        directory("results/binning/dastool/{sample}_bins"),
    threads: config["resources"]["threads_binning"]
    conda: "../envs/binning.yaml"
    log: "logs/dastool/{sample}.log"
    shell:
        """
        # Generate scaffold-to-bin tables
        Fasta_to_Contig2Bin.sh -i {input.metabat2} -e fa \
            > results/binning/dastool/{wildcards.sample}_metabat2.tsv
        Fasta_to_Contig2Bin.sh -i {input.concoct} -e fa \
            > results/binning/dastool/{wildcards.sample}_concoct.tsv
        Fasta_to_Contig2Bin.sh -i {input.maxbin2} -e fasta \
            > results/binning/dastool/{wildcards.sample}_maxbin2.tsv

        DAS_Tool \
            -i results/binning/dastool/{wildcards.sample}_metabat2.tsv,\
results/binning/dastool/{wildcards.sample}_concoct.tsv,\
results/binning/dastool/{wildcards.sample}_maxbin2.tsv \
            -l MetaBAT2,CONCOCT,MaxBin2 \
            -c {input.contigs} \
            -o results/binning/dastool/{wildcards.sample} \
            --write_bins \
            -t {threads} \
            --score_threshold 0.5 \
            --search_engine diamond \
            2> {log}

        mkdir -p {output}
        mv results/binning/dastool/{wildcards.sample}_DASTool_bins/* {output}/ \
            2>/dev/null || true
        """

rule collect_all_bins:
    """Collect refined bins from all samples"""
    input:
        expand("results/binning/dastool/{sample}_bins", sample=SAMPLES),
    output:
        directory("results/binning/dastool/all_bins"),
    shell:
        """
        mkdir -p {output}
        for d in {input}; do
            for f in "$d"/*.fa; do
                [ -f "$f" ] && cp "$f" {output}/
            done
        done
        """
