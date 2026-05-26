# =============================================================================
# Genome Binning Rules
# MetaBAT2, CONCOCT, MaxBin2 + DAS Tool integration
# =============================================================================

rule map_reads_to_assembly:
    """Map reads back to assembly for coverage estimation"""
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        r1="results/qc/{sample}_R1.clean.fastq.gz",
        r2="results/qc/{sample}_R2.clean.fastq.gz",
    output:
        bam="results/binning/mapping/{sample}.sorted.bam",
        bai="results/binning/mapping/{sample}.sorted.bam.bai",
    threads: 8
    log:
        "logs/binning/mapping_{sample}.log"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.contigs} {input.r1} {input.r2} 2> {log} \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule metabat2_depth:
    """Calculate contig depth for MetaBAT2"""
    input:
        bam="results/binning/mapping/{sample}.sorted.bam",
    output:
        depth="results/binning/metabat2/{sample}_depth.txt",
    log:
        "logs/binning/metabat2_depth_{sample}.log"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} \
            2> {log}
        """

rule metabat2_bin:
    """Genome binning with MetaBAT2"""
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        depth="results/binning/metabat2/{sample}_depth.txt",
    output:
        directory("results/binning/metabat2/{sample}_bins"),
    params:
        min_contig=config["binning"]["metabat2"]["min_contig"],
        max_p=config["binning"]["metabat2"]["max_p"],
        min_s=config["binning"]["metabat2"]["min_s"],
    threads: 4
    log:
        "logs/binning/metabat2_{sample}.log"
    shell:
        """
        mkdir -p {output}
        metabat2 \
            -i {input.contigs} \
            -a {input.depth} \
            -o {output}/bin \
            --minContig {params.min_contig} \
            --maxP {params.max_p} \
            --minS {params.min_s} \
            -t {threads} \
            2> {log}
        """

rule concoct_bin:
    """Genome binning with CONCOCT"""
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        bam="results/binning/mapping/{sample}.sorted.bam",
    output:
        directory("results/binning/concoct/{sample}_bins"),
    params:
        chunk_size=config["binning"]["concoct"]["chunk_size"],
    threads: 4
    log:
        "logs/binning/concoct_{sample}.log"
    shell:
        """
        mkdir -p results/binning/concoct/{wildcards.sample}_work
        cut_up_fasta.py {input.contigs} \
            -c {params.chunk_size} -o 0 \
            --merge_last \
            -b results/binning/concoct/{wildcards.sample}_work/contigs_10K.bed \
            > results/binning/concoct/{wildcards.sample}_work/contigs_10K.fa

        concoct_coverage_table.py \
            results/binning/concoct/{wildcards.sample}_work/contigs_10K.bed \
            {input.bam} \
            > results/binning/concoct/{wildcards.sample}_work/coverage_table.tsv

        concoct \
            --composition_file results/binning/concoct/{wildcards.sample}_work/contigs_10K.fa \
            --coverage_file results/binning/concoct/{wildcards.sample}_work/coverage_table.tsv \
            -b results/binning/concoct/{wildcards.sample}_work/ \
            -t {threads} \
            2> {log}

        merge_cutup_clustering.py \
            results/binning/concoct/{wildcards.sample}_work/clustering_gt1000.csv \
            > results/binning/concoct/{wildcards.sample}_work/clustering_merged.csv

        mkdir -p {output}
        extract_fasta_bins.py \
            {input.contigs} \
            results/binning/concoct/{wildcards.sample}_work/clustering_merged.csv \
            --output_path {output}
        """

rule maxbin2_bin:
    """Genome binning with MaxBin2"""
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        depth="results/binning/metabat2/{sample}_depth.txt",
    output:
        directory("results/binning/maxbin2/{sample}_bins"),
    params:
        min_contig=config["binning"]["maxbin2"]["min_contig"],
        prob=config["binning"]["maxbin2"]["prob_threshold"],
    threads: 4
    log:
        "logs/binning/maxbin2_{sample}.log"
    shell:
        """
        mkdir -p {output}
        # Extract abundance from MetaBAT2 depth file
        cut -f1,3 {input.depth} | tail -n +2 > results/binning/maxbin2/{wildcards.sample}_abund.txt

        run_MaxBin.pl \
            -contig {input.contigs} \
            -abund results/binning/maxbin2/{wildcards.sample}_abund.txt \
            -out {output}/bin \
            -min_contig_length {params.min_contig} \
            -prob_threshold {params.prob} \
            -thread {threads} \
            2> {log}
        """

rule dastool_refine:
    """Integrate binning results with DAS Tool"""
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        metabat2="results/binning/metabat2/{sample}_bins",
        concoct="results/binning/concoct/{sample}_bins",
        maxbin2="results/binning/maxbin2/{sample}_bins",
    output:
        directory("results/binning/dastool/{sample}_DASTool_bins"),
    params:
        score=config["binning"]["dastool"]["score_threshold"],
        engine=config["binning"]["dastool"]["search_engine"],
        prefix="results/binning/dastool/{sample}",
    threads: 8
    log:
        "logs/binning/dastool_{sample}.log"
    shell:
        """
        # Generate scaffold-to-bin files
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
            -l metabat2,concoct,maxbin2 \
            -c {input.contigs} \
            -o {params.prefix} \
            --score_threshold {params.score} \
            --search_engine {params.engine} \
            --threads {threads} \
            --write_bins \
            2> {log}
        """
