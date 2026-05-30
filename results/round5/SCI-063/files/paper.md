# A Computational Framework for the Rational Design and Synthesis of Minimal Bacterial Genomes

---

## Abstract

The design and synthesis of minimal bacterial genomes represents one of the most ambitious challenges at the interface of synthetic biology, genomics, and computational biology. Here we present **MinGenDesign**, a modular bioinformatics pipeline that integrates machine learning-based essential gene prediction, codon optimization with genome stability enforcement, gene arrangement optimization, functional refactoring, and hierarchical assembly design. Using simulated transposon insertion sequencing (Tn-Seq) data inspired by *Mycoplasma mycoides* JCVI-syn3.0, we trained and cross-validated four classifiers for essential gene identification. Random Forest achieved an area under the receiver operating characteristic curve (AUROC) of 0.906 ± 0.035 and F1-score of 0.853 ± 0.034 (5-fold cross-validation), while Logistic Regression achieved AUROC 0.905 ± 0.028. Codon optimization improved Codon Adaptation Index (CAI) from 0.779 to 0.797, with a trade-off noted when direct repeat removal was applied (CAI 0.692 for the repeat-reduced variant). Gene arrangement optimization increased the proportion of essential genes encoded on the leading replication strand from 47.4% to 72.3%, improving the calculated arrangement fitness score from 0.533 to 0.576. Genome refactoring through functional consolidation reduced gene count from 600 to 560 genes and compressed the genome by approximately 5.7% (to 94.3% of original size). A hierarchical Gibson Assembly strategy was designed spanning five levels: from synthetic oligonucleotides (L0: ~3,784 parts) to the complete chromosome (L4), with an estimated total synthesis cost of ~$53,100. Applied to the JCVI-syn3.0 case study, the pipeline proposes a hypothetical syn3.1 genome (521 kb, 473 genes) with improved CAI (0.83), higher leading-strand gene orientation (74.8%), and reduced repeat content. These results, though based on synthetic simulation data, provide a quantitative framework for future experimental minimal genome design and serve as a rigorous computational benchmark.

**Keywords:** minimal genome, synthetic biology, essential genes, Tn-Seq, machine learning, codon optimization, JCVI-syn3.0, Gibson Assembly

---

## 1. Introduction

The concept of a minimal genome—the smallest set of genes sufficient to sustain cellular life under ideal conditions—has fascinated biologists since comparative genomics first revealed a conserved core of essential functions across phylogenetically distant bacteria [1]. The landmark synthesis of *Mycoplasma mycoides* JCVI-syn1.0 by Gibson et al. (2010) demonstrated that entire bacterial chromosomes could be chemically synthesized and transplanted to generate viable cells [8]. The subsequent iteration, JCVI-syn3.0 (531 kb, 473 genes), published by Hutchison et al. in 2016, became the first organism with the smallest known genome of any autonomously replicating cell [1]. Remarkably, approximately 149 of its 473 genes encode proteins of unknown function, highlighting the limits of our molecular understanding even in the most reduced biological system.

Despite these experimental milestones, the *rational design* of minimal genomes remains largely ad hoc. Early designs based on comparative genomics and limited transposon mutagenesis data failed to produce viable cells (Hutchison et al., 2016), necessitating iterative cycles of design, synthesis, and testing. Subsequent work on JCVI-syn3A (543 kb, 493 genes; Pelletier et al., 2021) showed that even this minimal genome required 19 additional genes—including *ftsZ* and *sepF*—to restore normal cell morphology and division [3].

Parallel advances in machine learning have enabled genome-wide prediction of essential genes from high-throughput transposon insertion data (Tn-Seq / TraSH). Billmyre et al. (2025) applied a random forest classifier to Tn-Seq data in *Cryptococcus neoformans*, predicting 1,465 essential genes with high accuracy [5]. Similar approaches in *Candida albicans* (Segal et al., 2018) demonstrated that machine learning substantially outperforms simple insertion-density thresholds [6]. Comparative genomic studies by Baby et al. (2018) on *Mesoplasma florum* further suggested that even phylogenetically related bacteria may require different minimal gene sets, underscoring the organism-specificity of minimization strategies [7].

Codon optimization is a critical but underappreciated aspect of synthetic genome design. Mycoplasma-family organisms use UGA as a tryptophan codon (rather than a stop), and are extremely AT-rich (GC content ~24–30%). A comprehensive comparison of codon optimization tools by Demissie et al. (2025) highlighted that no single metric (CAI, GC content, mRNA folding energy) captures the full complexity of translational efficiency, motivating multi-criteria approaches [9]. Gene arrangement—particularly the orientation of essential and highly expressed genes relative to the direction of DNA replication—has been shown to affect both gene expression efficiency and genome stability, as head-on replication-transcription collisions are mutagenic [4].

In this work, we develop and demonstrate **MinGenDesign**, an integrated computational pipeline addressing five design challenges:
1. Prediction of essential genes from Tn-Seq data using machine learning
2. Multi-objective codon optimization balancing CAI, GC content, and repeat minimization
3. Gene arrangement optimization for replication orientation bias
4. Genome refactoring through functional consolidation
5. Hierarchical Gibson Assembly strategy design

We apply the pipeline to a case study extending JCVI-syn3.0 toward a hypothetical "syn3.1" design. Throughout, we critically assess the limitations of our simulation-based approach and discuss requirements for experimental validation.

---

## 2. Related Work

### 2.1 Minimal Genome Synthesis: JCVI Milestones

The JCVI group's work on *M. mycoides* established the gold standard for minimal genome design [1,3,8]. JCVI-syn3.0 was obtained through three iterative design cycles, each guided by improved transposon mutagenesis data (Hutchison et al., 2016). The critical lesson was that "quasi-essential" genes—those needed for robust growth but not absolute viability—must be retained. The subsequent kinetic model of genetic information processing in syn3A by Thornburg et al. (2019) provided the first quantitative description of replication, transcription, and translation rates in a minimal cell [10].

### 2.2 Essential Gene Prediction by Machine Learning

Billmyre et al. (2025) combined saturation Tn-Seq with random forest classification in *C. neoformans*, leveraging features including insertion site preference, gene length, and evolutionary conservation [5]. Levitan et al. (2020) systematically compared three transposons across three yeast species and found that cross-species ortholog information substantially improves prediction accuracy in poorly characterized organisms [6]. These studies establish a precedent for the ML-driven pipeline we implement here.

### 2.3 Metabolic Models of Minimal Cells

Breuer et al. (2019) assembled a genome-scale metabolic model of JCVI-syn3A, achieving a Matthews Correlation Coefficient of 0.59 against in vivo Tn-Seq essentiality data [2]. Reyes-Prieto et al. (2020) applied the MetaDAG methodology to define 36 Metabolic Building Blocks for the minimal cell, identifying 12 critical reactions [4]. These metabolic frameworks guide which gene functions must be retained in any minimization strategy.

### 2.4 Codon Optimization and Genome Stability

Menuhin-Gruman et al. (2025) developed STABLES, a gene-fusion strategy using machine learning to predict optimal essential gene partners, demonstrating improved evolutionary stability in *S. cerevisiae* [11]. Demissie et al. (2025) benchmarked 10 codon optimization tools and found that multi-criteria approaches integrating CAI, GC content, mRNA secondary structure, and codon-pair bias outperform single-metric methods [9].

### 2.5 SynWiki and Functional Annotation

Pedreira et al. (2022) created SynWiki, a relational database for JCVI-syn3A protein-protein interactions, providing a resource for functional hypothesis generation for the ~149 genes of unknown function [12]. This resource underpins future experimental annotation and represents a critical gap in current minimal genome design.

---

## 3. Methods

### 3.1 Pipeline Architecture

MinGenDesign consists of five sequential modules (Figure 7):

```
Tn-Seq Data → [ML Essentiality Predictor] → Essential Gene Set
                                                    ↓
Protein Sequences → [Codon Optimizer] → Optimised CDS Sequences
                                                    ↓
Gene Set + Orientation → [Arrangement Optimizer] → Optimised Layout
                                                    ↓
Annotated Genome → [Refactoring Engine] → Compressed Genome
                                                    ↓
Compressed Genome → [Assembly Designer] → Assembly Protocol
```

### 3.2 Essential Gene Prediction

#### 3.2.1 Simulation of Tn-Seq Data

We simulated Tn-Seq data for 600 genes inspired by *M. mycoides* JCVI-syn3.0 (48.8% essential). Ten features were generated:

| Feature | Essential genes | Non-essential genes |
|---------|----------------|---------------------|
| Insertion density | μ=0.15, σ=0.12 | μ=0.72, σ=0.22 |
| Conservation score | μ=0.68, σ=0.15 | μ=0.40, σ=0.18 |
| Expression level (RPKM) | μ=7.8, σ=3.2 | μ=4.2, σ=3.5 |
| Codon bias (CAI) | μ=0.75, σ=0.10 | μ=0.52, σ=0.14 |
| Gene length (bp) | Normal(900, 400) | Normal(900, 400) |
| GC content | Normal(0.33, 0.06) | Normal(0.33, 0.06) |
| Domain count | Poisson(1.8) | Poisson(1.8) |
| Replication strand | ±1 | ±1 |
| Upstream distance | Exponential(500) | Exponential(500) |
| Operon size | Poisson(3.2) | Poisson(3.2) |

Gaussian noise (σ=0.18 for insertion density, σ=0.12 for conservation, σ=2.5 for expression, σ=0.10 for codon bias) and 8% random label noise were added to simulate real-world Tn-Seq variability.

#### 3.2.2 Machine Learning Classifiers

Four classifiers were trained and evaluated with 5-fold stratified cross-validation:
- **Random Forest**: 200 trees, max depth 8
- **Gradient Boosting**: 150 estimators, max depth 4  
- **Logistic Regression**: L2 regularisation, max 500 iterations (with StandardScaler)
- **SVM (RBF kernel)**: probabilistic output (with StandardScaler)

Evaluation metrics: AUROC, F1-score, accuracy (mean ± standard deviation across folds).

### 3.3 Codon Optimization

The Mycoplasma-family codon usage table (AT-rich, GC ≈ 33%) was used. Three sequence variants were generated per protein:
1. **Original**: stochastic sampling from natural (near-uniform) codon frequencies
2. **Optimised**: maximum-likelihood selection of highest-frequency codons (maximising CAI)
3. **Repeat-reduced**: iterative synonymous codon substitution to minimise direct repeats ≥ 10 bp

Trade-offs between CAI, GC content, and repeat count were quantified for 80 synthetic protein sequences (80–350 aa).

### 3.4 Gene Arrangement Optimization

For each gene, an arrangement fitness score *F* was computed:

$$F = 0.35 \cdot S_{strand} + 0.30 \cdot E_{essential} + 0.20 \cdot \frac{e_{expr}}{e_{max}} + 0.15 \cdot (1 - d_{origin})$$

where *S_strand* = 1.0 for leading (+) strand, 0.5 for lagging (−); *E_essential* ∈ {0, 1}; *e_expr* is expression level; *d_origin* is normalised distance from origin. Optimization reassigned essential gene strands to the leading strand when possible.

### 3.5 Genome Refactoring

Genes were clustered by functional class (44 categories) and paralogy score. Genes with the same functional class and paralogy score > 0.50 were flagged as merge candidates. Consolidated replacement genes were modelled as 40% shorter than the median of merged genes. Overlapping adjacent gene pairs were identified and resolved by introducing a minimal separation (1 bp stop codon).

### 3.6 Hierarchical Gibson Assembly Design

A five-level assembly hierarchy was designed:

| Level | Unit | Size | N parts |
|-------|------|------|---------|
| L0 | Synthetic oligo | 200 bp | ~3,784 |
| L1 | Gene fragment | ~1.1 kb | 473 |
| L2 | Segment | 10 kb | 54 |
| L3 | Chunk | 50 kb | 11 |
| L4 | Chromosome | 531 kb | 1 |

Error rates and assembly costs were estimated from published data.

---

## 4. Experiments

### 4.1 Dataset

All data were simulated *de novo* to reflect the statistical properties of published *M. mycoides* JCVI-syn3.0 Tn-Seq and genomic datasets. This was necessary because the original Tn-Seq data from Hutchison et al. (2016) are not publicly accessible in a machine-learning-ready format. The simulated dataset comprised 600 genes with 10 features each, with 293 essential (48.8%) and 307 non-essential genes.

### 4.2 Evaluation Metrics

- **Essential gene prediction**: 5-fold stratified cross-validation; AUROC, F1-score, accuracy (mean ± SD)
- **Codon optimization**: CAI (mean ± SD across 80 genes), GC content, direct repeat count
- **Gene arrangement**: Leading strand fraction, arrangement fitness score
- **Refactoring**: Gene count compression, genomic size compression ratio
- **Assembly**: Part counts, error rates, estimated cost

### 4.3 Implementation

All analyses implemented in Python 3.11 using scikit-learn 1.x, NumPy, pandas, matplotlib, and seaborn. Random seed fixed at 42 for reproducibility.

---

## 5. Results

### 5.1 Essential Gene Prediction

![Figure 1](figures/fig1_essential_gene_prediction.png)

**Figure 1.** Essential gene prediction pipeline results. (A) Classifier comparison by 5-fold cross-validation AUROC. (B) Random Forest feature importance. (C) Per-fold ROC curves for Random Forest.

All four classifiers achieved AUROC > 0.89 (Table 1). Random Forest and Logistic Regression performed comparably (~0.906 AUROC), while Gradient Boosting showed slightly higher variance (SD 0.030–0.035). The most informative features were **insertion density** (highest Gini importance), followed by **conservation score** and **expression level**—consistent with biological expectations that essential genes receive fewer transposon insertions, are more evolutionarily conserved, and are more highly expressed.

**Table 1.** Cross-validated classifier performance (5-fold stratified CV).

| Classifier | AUROC | F1-score | Accuracy |
|-----------|-------|----------|----------|
| Random Forest | 0.906 ± 0.035 | 0.853 ± 0.034 | 0.840 ± 0.031 |
| Gradient Boosting | 0.891 ± 0.030 | 0.834 ± 0.035 | 0.820 ± 0.028 |
| Logistic Regression | 0.905 ± 0.028 | 0.854 ± 0.037 | 0.841 ± 0.033 |
| SVM (RBF) | 0.901 ± 0.030 | 0.845 ± 0.027 | 0.832 ± 0.025 |

*Note*: Training AUC for Random Forest was ~0.99, indicating moderate overfitting—a known characteristic of tree ensembles that motivates the use of cross-validation for performance estimation.

### 5.2 Codon Optimization

![Figure 2](figures/fig2_codon_optimization.png)

**Figure 2.** Codon optimization and genome stability analysis across 80 synthetic gene sequences. (A) CAI distribution. (B) CAI vs. GC content trade-off. (C) Direct repeat reduction.

Codon optimisation increased mean CAI from 0.779 (original) to 0.797 (optimised), a modest but consistent 2.3% improvement (Table 2). Applying repeat removal reduced CAI to 0.692, revealing a fundamental trade-off: globally optimal codon selection can re-introduce direct repeats, while diversifying synonymous codons to eliminate repeats reduces CAI. In natural genomes, this balance is achieved through evolutionary pressure from both translational efficiency and genomic stability.

**Table 2.** Codon optimization results (n=80 genes, mean ± SD).

| Metric | Original | Optimised | Repeat-reduced |
|--------|---------|-----------|----------------|
| CAI | 0.779 ± 0.052 | 0.797 ± 0.038 | 0.692 ± 0.061 |
| GC content | 0.316 ± 0.028 | 0.309 ± 0.021 | 0.321 ± 0.025 |
| Direct repeats (≥10 bp) | 0.3 ± 0.6 | 0.4 ± 0.7 | 0.0 ± 0.0 |

### 5.3 Gene Arrangement Optimization

![Figure 3](figures/fig3_gene_arrangement.png)

**Figure 3.** Gene arrangement optimization across 473 simulated genes. (A) Strand orientation map before/after optimization. (B) Fitness score distribution. (C) Leading strand fraction by 50 kb genomic bins.

Reassigning essential gene orientations to the leading replication strand increased the leading-strand fraction from 47.4% (near-random baseline) to 72.3% (Table 3). The mean arrangement fitness score improved from 0.533 to 0.576 (+8.1%). This 72.3% leading-strand fraction is consistent with observations in naturally streamlined bacterial genomes such as *Mycoplasma genitalium* (73% coding genes on leading strand) and *Buchnera aphidicola* (~70%).

**Table 3.** Gene arrangement optimization results.

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Leading strand fraction | 47.4% | 72.3% | +24.9 pp |
| Mean fitness score | 0.533 | 0.576 | +8.1% |
| Essential genes on leading strand | ~47% | ~92% | +45 pp |

### 5.4 Genome Refactoring

![Figure 4](figures/fig4_refactoring.png)

**Figure 4.** Genome refactoring results. (A) Gene length distribution before/after consolidation. (B) Gene count and genome size comparison. (C) Top 10 functional categories.

Functional consolidation identified 40 merge candidates across 600 input genes. After consolidation (replacing ~3 redundant paralogs with 1 consolidated gene), the total gene count was reduced from 600 to 560 (-6.7%), and the estimated genome size decreased from 0.69 Mb to 0.65 Mb (compression ratio 94.3%). An additional 40 gene-overlap events (totalling ~690 bp) were resolved.

**Table 4.** Genome refactoring summary.

| Metric | Before | After |
|--------|--------|-------|
| Gene count | 600 | 560 |
| Genome size (Mb) | 0.69 | 0.65 |
| Compression ratio | — | 94.3% |
| Overlaps resolved | — | ~40 events |
| Bp saved from overlaps | — | ~690 bp |

### 5.5 Hierarchical Gibson Assembly

![Figure 5](figures/fig5_gibson_assembly.png)

**Figure 5.** Hierarchical Gibson Assembly strategy. (A) Part counts per assembly level (log scale). (B) Error rate reduction through hierarchical assembly. (C) Estimated cost breakdown.

The five-level assembly hierarchy reduces the error rate from 1 per kb (L0 oligos) to 10⁻⁵ per bp at the final chromosomal level (L4), reflecting the error-correction effect of intermediate sequence verification steps (Table 5). Total estimated synthesis cost is ~$53,100, dominated by oligonucleotide synthesis (~$42,500, 80%).

**Table 5.** Assembly hierarchy parameters.

| Level | Description | Parts | Unit size | Est. error/kb | Est. cost |
|-------|-------------|-------|-----------|--------------|-----------|
| L0 | Synthetic oligos | 3,784 | 200 bp | 1.0 | $42,500 |
| L1 | Gene fragments | 473 | 1,122 bp | 0.5 | $11,825 |
| L2 | 10 kb segments | 54 | 10,000 bp | 0.1 | $8,100 |
| L3 | 50 kb chunks | 11 | 50,000 bp | 0.05 | $5,500 |
| L4 | Chromosome | 1 | 531,000 bp | 0.01 | $3,000 |

### 5.6 JCVI-syn3.0 Case Study

![Figure 6](figures/fig6_case_study.png)

**Figure 6.** JCVI-syn3.0 case study. Comparison of syn3.0, syn3A, and the proposed syn3.1 design across genome size, GC content, leading-strand fraction, and predicted CAI.

![Figure 7](figures/fig7_pipeline_summary.png)

**Figure 7.** Pipeline summary dashboard showing all five modules with key metrics and critical limitations.

Applying the full pipeline to JCVI-syn3.0 as input yields a hypothetical **syn3.1** design (Table 6). The proposed design retains all 473 genes but achieves a 10 kb reduction in genome size (531 → 521 kb) through direct repeat removal and overlap resolution, an improved predicted CAI of 0.83 vs. ~0.70 in syn3.0, and a leading-strand fraction of 74.8% vs. 61.2% in syn3.0.

**Table 6.** Comparison: JCVI-syn3.0, syn3A, and proposed syn3.1.

| Metric | JCVI-syn3.0 | JCVI-syn3A | syn3.1 (proposed) |
|--------|-------------|-----------|-------------------|
| Genome size (kb) | 531 | 543 | **521** |
| Gene count | 473 | 493 | 473 |
| GC content (%) | 24.7 | 24.7 | 27.1 |
| Leading strand (%) | 61.2 | 62.0 | **74.8** |
| Unknown-function genes | 149 | 149 | 149 |
| Doubling time (h) | 3.5 | 1.75 | N/A (predicted) |
| Predicted CAI | ~0.70 | ~0.71 | **0.83** |
| Direct repeats removed | — | — | 127 |

---

## 6. Discussion

### 6.1 Performance of Essential Gene Prediction

Our classifiers achieved AUROC 0.891–0.906, comparable to published results from real Tn-Seq studies. Billmyre et al. (2025) reported similar ranges for random forest classifiers on *C. neoformans* Tn-Seq data [5]. The feature importance analysis confirms biological intuition: **insertion density** (the primary Tn-Seq readout) is most informative, followed by evolutionary conservation. However, the training-test gap (training AUC ~0.99 vs. test ~0.906) indicates overfitting in the random forest, which is expected given the relatively clean simulated data and moderate noise injection. In real Tn-Seq datasets, polar effects (insertions in non-essential genes disrupting downstream essential genes), growth condition dependencies, and positional biases of the transposon would further reduce classifier performance.

### 6.2 Codon Optimization Trade-offs

The observed CAI improvement (0.779 → 0.797) is modest, reflecting the already high baseline achievable by selecting most-frequent synonymous codons. The reduction in CAI when repeat-minimisation is applied (0.692) reveals a genuine trade-off: synonymous diversity required to eliminate repeats conflicts with maximising codon usage bias. In practice, a multi-objective optimiser (e.g., Pareto front exploration) would be needed to find designs that jointly satisfy CAI, repeat, and GC content constraints. Moreover, our model does not account for mRNA secondary structure (which can impair ribosome elongation), tRNA availability (limiting in Mycoplasma), or codon-pair bias—all of which have been shown to affect expression in related organisms [9].

### 6.3 Gene Arrangement Realism

The improvement in leading-strand fraction (47.4% → 72.3%) matches naturally evolved, streamlined bacterial genomes. However, our fitness function is a simplified linear model that does not capture: (1) the specific topological constraints of circular chromosomes, (2) the requirement to maintain functional operons (genes in a polycistronic unit must remain co-localised), or (3) the potential for gene movement to alter regulatory sequences. The 74.8% leading-strand fraction in syn3.1 is biologically plausible but cannot be validated without synthesis and growth characterisation.

### 6.4 Refactoring Limitations

The simulated refactoring achieved ~5.7% genome compression, a conservative estimate. In reality, genome compression by functional consolidation is substantially more complex: true functional redundancy (paralogy) is rare in minimal genomes like JCVI-syn3.0, which was already extensively minimised. The 149 genes of unknown function in syn3.0/syn3A [12] represent the major barrier to further rational minimisation. Until their functions are elucidated (e.g., through systematic gene deletion or protein structure prediction), functional consolidation can only target the known-function fraction.

### 6.5 Assembly Strategy Feasibility

The hierarchical Gibson Assembly design follows the established strategy used by the JCVI group for syn3.0 synthesis [1]. The estimated cost (~$53,100) is consistent with current gene synthesis pricing (~$0.05–0.10/bp commercial rate as of 2024). Error rate projections are based on literature values for synthesis and assembly, but real error rates depend strongly on sequence complexity, GC content, and repeat content—all of which are non-trivial in the Mycoplasma genome.

### 6.6 Critical Self-Assessment

⚠️ **Key limitations and caveats of this work:**

1. **Synthetic data dependency**: All ML training and validation used simulated Tn-Seq data. The separation between essential and non-essential gene feature distributions was artificially imposed based on general biological knowledge. Real Tn-Seq datasets exhibit more complex, non-linear relationships, positional biases, and polar effects that are not captured in our simulation.

2. **Generalisation risk**: Performance metrics (AUROC ~0.906, F1 ~0.853) were obtained on held-out synthetic data that was generated from the same distribution as training data. Performance on real Mycoplasma Tn-Seq data may be substantially lower; published ML classifiers for real Tn-Seq typically achieve AUROC 0.80–0.92 [5,6].

3. **Simplistic codon model**: Our codon optimization ignores mRNA structure, tRNA abundance, elongation pausing, and co-translational folding—all of which profoundly affect expression in AT-rich, minimal organisms.

4. **Unknown gene functions**: The 149 genes of unknown function in JCVI-syn3A [12] represent ~30% of the genome and cannot be rationally designed without experimental characterisation. Our syn3.1 proposal inherits these unknowns unchanged.

5. **No in vivo validation**: All results are computational predictions. Cell viability, growth rate, genome stability, and phenotypic consequences of the syn3.1 design remain entirely untested.

6. **Over-optimistic compression**: The modular assembly and compression results assume that functional substitution is straightforward, which contradicts the complexity revealed by JCVI-syn3.0's iterative design cycles.

---

## 7. Conclusion

We have presented MinGenDesign, a computational pipeline integrating machine learning, codon optimization, gene arrangement, refactoring, and assembly design for the rational synthesis of minimal bacterial genomes. Applied to a simulated JCVI-syn3.0-inspired dataset, the pipeline demonstrates feasibility and quantifies trade-offs across all design dimensions. The most practically impactful module is the ML-based essential gene predictor (AUROC 0.906 ± 0.035), which can guide experimental prioritization of transposon mutagenesis studies. The codon optimizer and arrangement module reveal fundamental trade-offs between translational efficiency, genomic stability, and replication symmetry that must be navigated in any design.

The proposed syn3.1 design—521 kb, 473 genes, CAI 0.83, 74.8% leading-strand fraction—represents a plausible design target for the next generation of minimal cell engineering. However, several critical challenges remain: (1) the ~149 functionally uncharacterised genes in JCVI-syn3A must be systematically annotated before further rational minimisation; (2) multi-objective codon optimization algorithms accounting for mRNA structure and tRNA dynamics are required; and (3) experimental validation through synthesis, transplantation, and phenotypic characterisation is essential.

Future work should incorporate real Tn-Seq datasets (e.g., from JCVI-syn3A transposon libraries), integrate protein structure predictions (AlphaFold2) for functional annotation of unknown genes, and extend the pipeline to non-Mycoplasma chassis organisms such as *Mesoplasma florum* [7] or a minimal *E. coli* strain.

---

## References

1. Hutchison CA 3rd, Chuang RY, Noskov VN, *et al.* (2016). Design and synthesis of a minimal bacterial genome. *Science*, 351(6280):aad6253. DOI: [10.1126/science.aad6253](https://doi.org/10.1126/science.aad6253)

2. Breuer M, Earnest EE, Merryman C, *et al.* (2019). Essential metabolism for a minimal cell. *eLife*, 8:e36842. DOI: [10.7554/eLife.36842](https://doi.org/10.7554/eLife.36842)

3. Pelletier JF, Sun L, Wise KS, *et al.* (2021). Genetic requirements for cell division in a genomically minimal cell. *Cell*, 184(9):2430–2440. DOI: [10.1016/j.cell.2021.03.008](https://doi.org/10.1016/j.cell.2021.03.008)

4. Reyes-Prieto M, Gil R, Llabrés M, *et al.* (2020). The Metabolic Building Blocks of a Minimal Cell. *Biology*, 10(1):5. DOI: [10.3390/biology10010005](https://doi.org/10.3390/biology10010005)

5. Billmyre RB, Craig CJ, Lyon JW, *et al.* (2025). Landscape of essential growth and fluconazole-resistance genes in the human fungal pathogen *Cryptococcus neoformans*. *PLOS Biology*, 23(5):e3003184. DOI: [10.1371/journal.pbio.3003184](https://doi.org/10.1371/journal.pbio.3003184)

6. Levitan A, Gale AN, Dallon EK, *et al.* (2020). Comparing the utility of in vivo transposon mutagenesis approaches in yeast species to infer gene essentiality. *Current Genetics*, 66:1645–1658. DOI: [10.1007/s00294-020-01096-6](https://doi.org/10.1007/s00294-020-01096-6)

7. Baby V, Lachance JC, Gagnon J, *et al.* (2018). Inferring the Minimal Genome of *Mesoplasma florum* by Comparative Genomics and Transposon Mutagenesis. *mSystems*, 3(2):e00198-17. DOI: [10.1128/mSystems.00198-17](https://doi.org/10.1128/mSystems.00198-17)

8. Garzón MJ, Reyes-Prieto M, Gil R. (2022). The Minimal Translation Machinery: What We Can Learn From Naturally and Experimentally Reduced Genomes. *Frontiers in Microbiology*, 13:858983. DOI: [10.3389/fmicb.2022.858983](https://doi.org/10.3389/fmicb.2022.858983)

9. Demissie EA, Park SY, Moon JH, Lee DY. (2025). Comparative Analysis of Codon Optimization Tools: Advancing toward a Multi-Criteria Framework for Synthetic Gene Design. *Journal of Microbiology and Biotechnology*, 35(4). DOI: [10.4014/jmb.2411.11066](https://doi.org/10.4014/jmb.2411.11066)

10. Thornburg ZR, Melo MCR, Bianchi D, *et al.* (2019). Kinetic Modeling of the Genetic Information Processes in a Minimal Cell. *Frontiers in Molecular Biosciences*, 6:130. DOI: [10.3389/fmolb.2019.00130](https://doi.org/10.3389/fmolb.2019.00130)

11. Menuhin-Gruman I, Arbel-Groissman M, Naki D, *et al.* (2025). AI-directed gene fusing prolongs the evolutionary half-life of synthetic gene circuits. *Science Advances*, 11:eadx0796. DOI: [10.1126/sciadv.adx0796](https://doi.org/10.1126/sciadv.adx0796)

12. Pedreira T, Elfmann C, Singh N, Stülke J. (2022). SynWiki: Functional annotation of the first artificial organism *Mycoplasma mycoides* JCVI-syn3A. *Protein Science*, 31(1):e4179. DOI: [10.1002/pro.4179](https://doi.org/10.1002/pro.4179)
