# ============================================================
# Step 6: Multivariate Statistical Analysis
#         Gut microbiota–disease association
# ============================================================

rule alpha_diversity:
    """Alpha diversity: Shannon, Simpson, observed features, Chao1"""
    input:
        profiles="results/taxonomy/merged_metaphlan4.tsv",
        metadata="config/samples.tsv",
    output:
        table="results/statistics/alpha_diversity.tsv",
        figure="figures/alpha_diversity_boxplot.svg",
    params:
        metrics=config["statistics"]["alpha_diversity_metrics"],
    conda: "../envs/statistics.yaml"
    log: "logs/alpha_diversity.log"
    script:
        "../scripts/alpha_diversity.py"

rule beta_diversity:
    """Beta diversity: PCoA on Bray-Curtis, Jaccard, Aitchison distances"""
    input:
        profiles="results/taxonomy/merged_metaphlan4.tsv",
        metadata="config/samples.tsv",
    output:
        pcoa="results/statistics/beta_diversity_pcoa.tsv",
        permanova="results/statistics/permanova_results.tsv",
        figure="figures/beta_diversity_pcoa.svg",
    params:
        metrics=config["statistics"]["beta_diversity_metrics"],
    conda: "../envs/statistics.yaml"
    log: "logs/beta_diversity.log"
    script:
        "../scripts/beta_diversity.py"

rule taxonomic_barplot:
    """Stacked barplot of taxonomic composition at phylum/genus level"""
    input:
        profiles="results/taxonomy/merged_metaphlan4.tsv",
        metadata="config/samples.tsv",
    output:
        "figures/taxonomic_barplot.svg",
    conda: "../envs/statistics.yaml"
    log: "logs/taxonomic_barplot.log"
    script:
        "../scripts/taxonomic_barplot.py"

rule differential_abundance:
    """Differential abundance: ALDEx2 + ANCOM-BC + MaAsLin 2 consensus"""
    input:
        profiles="results/taxonomy/merged_metaphlan4.tsv",
        pathways="results/functional/humann3/merged_pathabundance_cpm.tsv",
        metadata="config/samples.tsv",
    output:
        da_table="results/statistics/differential_abundance.tsv",
        maaslin2_dir=directory("results/statistics/maaslin2_results"),
    params:
        fdr=config["statistics"]["differential_abundance"]["fdr_threshold"],
        effect=config["statistics"]["differential_abundance"]["effect_size_threshold"],
        min_prev=config["statistics"]["min_prevalence"],
        min_abund=config["statistics"]["min_abundance"],
    conda: "../envs/statistics.yaml"
    log: "logs/differential_abundance.log"
    script:
        "../scripts/differential_abundance.py"

rule functional_heatmap:
    """Heatmap of differentially abundant pathways"""
    input:
        pathways="results/functional/humann3/merged_pathabundance_cpm.tsv",
        da="results/statistics/differential_abundance.tsv",
        metadata="config/samples.tsv",
    output:
        "figures/functional_heatmap.svg",
    conda: "../envs/statistics.yaml"
    log: "logs/functional_heatmap.log"
    script:
        "../scripts/functional_heatmap.py"

rule mag_quality_plot:
    """Scatter plot: MAG completeness vs contamination"""
    input:
        report="results/mag_quality/checkm2_report.tsv",
    output:
        "figures/mag_quality_scatter.svg",
    conda: "../envs/statistics.yaml"
    log: "logs/mag_quality_plot.log"
    script:
        "../scripts/mag_quality_plot.py"
