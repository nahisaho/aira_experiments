# =============================================================================
# Statistical Analysis Rules
# Alpha/Beta diversity, PERMANOVA, differential abundance
# =============================================================================

rule alpha_diversity:
    """Calculate alpha diversity metrics"""
    input:
        "results/taxonomy/merged_metaphlan_profiles.tsv",
    output:
        "results/stats/alpha_diversity.tsv",
    params:
        metrics=config["statistics"]["alpha_diversity_metrics"],
    log:
        "logs/stats/alpha_diversity.log"
    script:
        "../scripts/alpha_diversity.py"

rule beta_diversity:
    """Calculate beta diversity and PCoA"""
    input:
        "results/taxonomy/merged_metaphlan_profiles.tsv",
    output:
        pcoa="results/stats/beta_diversity_pcoa.tsv",
        dm="results/stats/distance_matrices/",
    params:
        metrics=config["statistics"]["beta_diversity_metrics"],
    log:
        "logs/stats/beta_diversity.log"
    script:
        "../scripts/beta_diversity.py"

rule permanova:
    """PERMANOVA test for group differences"""
    input:
        dm="results/stats/distance_matrices/",
        metadata="config/config.yaml",
    output:
        "results/stats/permanova_results.tsv",
    params:
        permutations=config["statistics"]["multivariate"]["permanova_permutations"],
    log:
        "logs/stats/permanova.log"
    script:
        "../scripts/permanova.py"

rule differential_abundance:
    """Differential abundance analysis with MaAsLin2-style approach"""
    input:
        profiles="results/taxonomy/merged_metaphlan_profiles.tsv",
        pathways="results/functional/merged_pathabundance.tsv",
        metadata="config/config.yaml",
    output:
        taxa="results/stats/differential_abundance.tsv",
        outdir=directory("results/stats/maaslin2_results/"),
    params:
        normalization=config["statistics"]["multivariate"]["maaslin2_normalization"],
        transform=config["statistics"]["multivariate"]["maaslin2_transform"],
        min_prevalence=config["statistics"]["multivariate"]["min_prevalence"],
        min_abundance=config["statistics"]["multivariate"]["min_abundance"],
    log:
        "logs/stats/differential_abundance.log"
    script:
        "../scripts/differential_abundance.py"

rule plot_taxonomy_barplot:
    """Generate taxonomy composition barplot"""
    input:
        "results/taxonomy/merged_metaphlan_profiles.tsv",
    output:
        "figures/taxonomy_barplot.png",
    log:
        "logs/stats/plot_taxonomy.log"
    script:
        "../scripts/plot_taxonomy.py"

rule plot_alpha_diversity:
    """Generate alpha diversity boxplots"""
    input:
        "results/stats/alpha_diversity.tsv",
    output:
        "figures/alpha_diversity_boxplot.png",
    log:
        "logs/stats/plot_alpha.log"
    script:
        "../scripts/plot_alpha_diversity.py"

rule plot_beta_diversity:
    """Generate beta diversity PCoA plot"""
    input:
        "results/stats/beta_diversity_pcoa.tsv",
    output:
        "figures/beta_diversity_pcoa.png",
    log:
        "logs/stats/plot_beta.log"
    script:
        "../scripts/plot_beta_diversity.py"

rule plot_functional_heatmap:
    """Generate functional pathway heatmap"""
    input:
        "results/functional/merged_pathabundance.tsv",
    output:
        "figures/functional_heatmap.png",
    log:
        "logs/stats/plot_functional.log"
    script:
        "../scripts/plot_functional_heatmap.py"

rule plot_mag_quality:
    """Generate MAG quality scatter plot"""
    input:
        "results/mags/all_mags_summary.tsv",
    output:
        "figures/mag_quality_scatter.png",
    log:
        "logs/stats/plot_mag.log"
    script:
        "../scripts/plot_mag_quality.py"

rule plot_volcano:
    """Generate differential abundance volcano plot"""
    input:
        "results/stats/differential_abundance.tsv",
    output:
        "figures/differential_abundance_volcano.png",
    log:
        "logs/stats/plot_volcano.log"
    script:
        "../scripts/plot_volcano.py"
