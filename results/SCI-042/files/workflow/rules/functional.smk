# ============================================================
# Step 3: Functional Annotation
#         HUMAnN 3 + eggNOG-mapper integration
# ============================================================

rule humann3_profile:
    """HUMAnN 3: gene family & pathway abundance profiling"""
    input:
        r1="results/qc/hostfree/{sample}_R1.fastq.gz",
        r2="results/qc/hostfree/{sample}_R2.fastq.gz",
        metaphlan="results/taxonomy/metaphlan4/{sample}_profile.tsv",
    output:
        genefamilies="results/functional/humann3/{sample}_genefamilies.tsv",
        pathabundance="results/functional/humann3/{sample}_pathabundance.tsv",
        pathcoverage="results/functional/humann3/{sample}_pathcoverage.tsv",
    params:
        nuc_db=config["databases"]["humann3_db"]["nucleotide"],
        prot_db=config["databases"]["humann3_db"]["protein"],
        outdir="results/functional/humann3",
    threads: config["resources"]["threads_annotation"]
    resources:
        mem_mb=64000,
    conda: "../envs/functional.yaml"
    log: "logs/humann3/{sample}.log"
    shell:
        """
        # Concatenate paired reads for HUMAnN 3
        cat {input.r1} {input.r2} > /tmp/{wildcards.sample}_concat.fastq.gz

        humann \
            --input /tmp/{wildcards.sample}_concat.fastq.gz \
            --output {params.outdir}/{wildcards.sample}_tmp \
            --taxonomic-profile {input.metaphlan} \
            --nucleotide-database {params.nuc_db} \
            --protein-database {params.prot_db} \
            --threads {threads} \
            --search-mode uniref90 \
            --memory-use maximum \
            2> {log}

        # Move and rename outputs
        mv {params.outdir}/{wildcards.sample}_tmp/*_genefamilies.tsv {output.genefamilies}
        mv {params.outdir}/{wildcards.sample}_tmp/*_pathabundance.tsv {output.pathabundance}
        mv {params.outdir}/{wildcards.sample}_tmp/*_pathcoverage.tsv {output.pathcoverage}
        rm -rf {params.outdir}/{wildcards.sample}_tmp /tmp/{wildcards.sample}_concat.fastq.gz
        """

rule humann3_normalize:
    """Normalize HUMAnN 3 output to CPM"""
    input:
        gf="results/functional/humann3/{sample}_genefamilies.tsv",
        pa="results/functional/humann3/{sample}_pathabundance.tsv",
    output:
        gf="results/functional/humann3/{sample}_genefamilies_cpm.tsv",
        pa="results/functional/humann3/{sample}_pathabundance_cpm.tsv",
    conda: "../envs/functional.yaml"
    log: "logs/humann3_norm/{sample}.log"
    shell:
        """
        humann_renorm_table -i {input.gf} -o {output.gf} \
            --units cpm --special n 2>> {log}
        humann_renorm_table -i {input.pa} -o {output.pa} \
            --units cpm --special n 2>> {log}
        """

rule humann3_merge:
    """Merge normalized HUMAnN 3 tables across samples"""
    input:
        gf=expand("results/functional/humann3/{sample}_genefamilies_cpm.tsv", sample=SAMPLES),
        pa=expand("results/functional/humann3/{sample}_pathabundance_cpm.tsv", sample=SAMPLES),
    output:
        gf="results/functional/humann3/merged_genefamilies_cpm.tsv",
        pa="results/functional/humann3/merged_pathabundance_cpm.tsv",
    conda: "../envs/functional.yaml"
    log: "logs/humann3_merge.log"
    shell:
        """
        humann_join_tables -i results/functional/humann3/ \
            --file_name genefamilies_cpm -o {output.gf} 2>> {log}
        humann_join_tables -i results/functional/humann3/ \
            --file_name pathabundance_cpm -o {output.pa} 2>> {log}
        """

rule humann3_regroup_ko:
    """Regroup gene families to KEGG Orthologs"""
    input:
        "results/functional/humann3/merged_genefamilies_cpm.tsv",
    output:
        "results/functional/humann3/merged_ko_cpm.tsv",
    conda: "../envs/functional.yaml"
    log: "logs/humann3_regroup.log"
    shell:
        """
        humann_regroup_table -i {input} -o {output} \
            -g uniref90_ko 2> {log}
        """

rule eggnog_mapper:
    """eggNOG-mapper: ORF prediction + functional annotation on assembled contigs"""
    input:
        contigs="results/assembly/{sample}/final_contigs.fa",
    output:
        "results/functional/eggnog/{sample}.emapper.annotations",
    params:
        db=config["databases"]["eggnog_db"],
        prefix="results/functional/eggnog/{sample}",
    threads: config["resources"]["threads_annotation"]
    conda: "../envs/functional.yaml"
    log: "logs/eggnog/{sample}.log"
    shell:
        """
        emapper.py \
            -i {input.contigs} \
            --itype metagenome \
            -m diamond \
            --data_dir {params.db} \
            --output {params.prefix} \
            --cpu {threads} \
            --override \
            --dbmem \
            --go_evidence all \
            --pfam_realign realign \
            2> {log}
        """
