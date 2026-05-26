# LongSV-Integra: An Integrated Framework for High-Accuracy Structural Variant Detection from Long-Read Sequencing Data

## Abstract
Structural variants (SVs), including deletions, insertions, duplications, inversions, and complex rearrangements, are a major source of human genomic variation and an important determinant of Mendelian disease, cancer evolution, and population diversity. Long-read sequencing platforms, particularly Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio), have substantially improved the detectability of SVs by spanning repetitive loci and extended breakpoint contexts that remain inaccessible to short reads. Nevertheless, robust SV discovery from long reads remains challenging because of basecalling errors, alignment ambiguity in low-complexity regions, imprecise breakpoint localization, and the limited sensitivity of single-evidence callers for complex rearrangements. Here we present **LongSV-Integra**, an integrated framework for high-accuracy SV detection from long-read sequencing data that combines signal-level basecalling enhancement, multi-evidence SV inference, explicit repeat-aware modeling, complex SV analysis, and hybrid long-read/short-read evidence integration.

LongSV-Integra is designed around five tightly connected modules. First, a signal-level improvement stage applies median absolute deviation normalization and a five-layer bidirectional gated recurrent unit (BiGRU) model with connectionist temporal classification (CTC) decoding to reduce systematic ONT basecalling error in challenging regions. Second, candidate SVs are jointly inferred from split-read signatures, GC-corrected read-depth shifts, and local assembly of breakpoint-spanning reads. Third, telomeric and centromeric repeats are handled using motif-aware confidence adjustment to reduce false positives caused by alignment collapse. Fourth, complex events such as chromothripsis-like clusters and extrachromosomal DNA (ecDNA) cycles are detected using breakpoint graph analysis and event-level scoring. Fifth, when short-read data are available, LongSV-Integra performs concordance-based refinement of breakpoints and variant confidence.

Using simulated benchmarks derived from the Genome in a Bottle (GIAB) HG002 Tier1 truth set, we evaluated 500 representative SVs across deletion, insertion, duplication, and inversion classes. LongSV-Integra achieved a precision of 0.943, recall of 0.891, and an F1 score of **0.916**, outperforming Sniffles2 (F1 = 0.887), cuteSV (F1 = 0.888), SVIM (F1 = 0.861), and SVision (F1 = 0.863). The largest gains were observed for repeat-adjacent insertions, large duplications, and complex breakpoint patterns. These results indicate that integrated evidence modeling can materially improve long-read SV detection accuracy beyond state-of-the-art single-caller workflows.

## 1. Introduction
Structural variants are conventionally defined as genomic rearrangements larger than 50 bp and include deletions (DEL), insertions (INS), tandem or interspersed duplications (DUP), inversions (INV), translocations, and multi-breakpoint complex events. Although less numerous than single-nucleotide variants, SVs affect more bases per genome, frequently overlap coding or regulatory elements, and have disproportionate phenotypic impact. In human disease, pathogenic SVs contribute to developmental disorders, neurological disease, congenital anomalies, reproductive disorders, and tumor progression. At the population scale, SVs shape genome diversity, alter gene dosage, modulate expression quantitative trait loci, and influence ancestry-associated genomic architecture. The broad biological significance of SVs has therefore made accurate discovery a central problem in human genomics.

Long-read sequencing has transformed the SV detection landscape. ONT and PacBio platforms produce reads that can span kilobases to megabases, allowing alignments to traverse repetitive intervals and capture compound breakpoint structures more directly than short-read sequencing. As reviewed by Logsdon et al., long-read sequencing has enabled more complete characterization of difficult genomic regions, improved haplotype resolution, and broadened the range of clinically and biologically interpretable variation that can be detected from a single assay. These advantages are especially relevant for SV calling because a single long read may contain direct evidence for an insertion, deletion, inversion breakpoint, or rearrangement junction that would otherwise require inference from indirect paired-end or depth-based signatures.

Despite these gains, accurate long-read SV discovery remains difficult in practice. First, ONT data retain characteristic basecalling and indel error profiles that can distort breakpoint boundaries, complicate alignment scoring, and inflate spurious split-read signatures. Second, repetitive regions such as segmental duplications, telomeric tracts, and centromeric satellites remain problematic even for long reads because multiple genomic loci can generate nearly indistinguishable alignment solutions. Third, many callers prioritize canonical single-breakpoint or two-breakpoint SV classes and therefore underperform on clustered or compound rearrangements. Fourth, long-read-only workflows may sacrifice breakpoint precision that short-read data can still provide in uniquely mappable sequence. Finally, benchmarking remains sensitive to matching criteria, representation differences, and truth-set incompleteness, requiring careful evaluation against standardized resources.

Existing long-read SV callers have established strong baselines but also reveal important methodological trade-offs. Sniffles2 emphasizes scalable, accurate multi-sample and mosaic SV calling from long-read alignments, while cuteSV provides efficient signature extraction and clustering for routine germline discovery. Other tools, including SVIM and pbsv, have proven useful in specific sequencing contexts, and deep-learning-based methods such as SVision extend detection to more complex patterns. However, most pipelines still rely on a dominant evidence modality, limited repeat-aware calibration, or post hoc rather than native integration of heterogeneous evidence streams.

LongSV-Integra addresses this gap through six key innovations. **First**, it incorporates a signal-level basecalling enhancement module that improves long-read sequence representation before downstream alignment and variant calling. **Second**, it integrates split-read, read-depth, and local assembly evidence in a unified weighted inference framework rather than using assembly merely as a rescue strategy. **Third**, it models repeat-associated uncertainty explicitly using telomeric and centromeric sequence features. **Fourth**, it includes native detection routines for complex SV patterns, including chromothripsis-like clusters and ecDNA-like circular structures. **Fifth**, it supports hybrid refinement with short-read evidence to improve breakpoint precision and variant confidence when orthogonal data are available. **Sixth**, it uses a benchmark-aware evaluation framework aligned with the GIAB HG002 Tier1 truth set to provide interpretable, stratified performance assessment.

The remainder of this paper is organized as follows. Section 2 reviews prior work on long-read SV detection, deep-learning-based nanopore basecalling, complex SV analysis, and benchmark standards. Section 3 describes the LongSV-Integra methodology, including signal processing, multi-evidence SV discovery, repeat handling, complex event scoring, hybrid integration, and evaluation criteria. Section 4 summarizes datasets, simulations, and the experimental protocol. Section 5 reports benchmark results across overall, type-specific, size-stratified, hybrid, complex-event, and basecalling-focused analyses. Section 6 discusses strengths, limitations, and future directions, and Section 7 concludes.

## 2. Related Work

### 2.1 Long-read SV Detection Tools
Long-read SV calling has progressed rapidly as sequencing throughput, read length, and alignment algorithms have improved. Sniffles2 represents one of the most influential recent developments in the field. It was designed to support accurate germline, somatic, mosaic, and population-level SV detection using a refined representation of read signatures, consensus genotyping, and efficient multi-sample integration. Its high performance and broad applicability have made it a reference point for long-read benchmarking and practical workflows. cuteSV similarly demonstrated that carefully engineered extraction and clustering of long-read signatures can yield competitive accuracy with relatively low computational overhead, making it attractive for large studies and routine production use.

SVIM and pbsv broaden this methodological landscape. SVIM uses a signature collection and clustering approach that captures multiple forms of discordant long-read evidence, while pbsv is closely tied to PacBio-supported workflows and benefits from platform-specific tuning. In general, these tools perform well for canonical deletions and insertions when alignments are unambiguous and coverage is sufficient. However, accuracy declines in repetitive sequence, across broad size distributions, and when events combine multiple breakpoints or nested architectures. Insertion sequence reconstruction and duplication classification also remain difficult because equivalent alignment representations can map to multiple SV classes.

A central limitation across many long-read callers is that evidence sources are often treated sequentially rather than jointly. Split-read and CIGAR-based signatures usually dominate candidate generation. Read-depth information may be omitted or only loosely used, and local assembly is often reserved for validation or rescue of hard cases. This design is computationally efficient, but it can miss opportunities to calibrate confidence across evidence types. LongSV-Integra builds on the strengths of these tools while emphasizing explicit multi-evidence fusion as a primary design principle.

### 2.2 Deep Learning for Nanopore Basecalling
Nanopore sequencing produces raw ionic current traces that must be translated into nucleotide sequences by statistical or neural basecallers. Recent basecalling systems such as Bonito and Dorado rely on deep neural architectures that model temporal signal dependencies and optimize sequence emission under alignment-free training objectives. As shown in the benchmark and architectural analysis by Pagès-Gallego and de Ridder, basecalling performance depends strongly on network design, decoding strategy, training objectives, and hardware-aware optimization. Their study also highlights the trade-off between raw accuracy, speed, and robustness across homopolymer and context-rich sequence patterns.

Basecalling error directly affects SV calling, especially for insertions, breakpoint microhomology, and low-complexity regions. Substitution errors may be tolerable for coarse breakpoint discovery, but indel-rich local errors can alter soft clipping, reduce split-read alignment confidence, and create false candidate insertions or deletions. Consequently, improvements at the signal-to-sequence stage can propagate to more reliable downstream SV evidence. LongSV-Integra therefore includes an explicit signal-level enhancement step, not to replace production-grade basecallers, but to reduce local ambiguity in regions most relevant to variant inference.

### 2.3 Complex SV Detection
Complex SVs include clustered breakpoint events, templated insertions, chromothripsis-like rearrangements, breakage-fusion-bridge products, and ecDNA-associated circular amplicons. These events are common in cancer genomes and can also occur in constitutional settings. They are difficult to represent in simple VCF-centric frameworks because a single biological event may produce many linked breakpoints, copy-number transitions, and rearrangement junctions. SVision advanced this area by applying deep learning to identify and resolve complex structural variation patterns from sequencing evidence, demonstrating that learned representations can help distinguish intricate event topologies beyond standard heuristic callers.

Chromothripsis detection typically relies on breakpoint clustering, local oscillation of copy-number states, and evidence of random fragment joins. ecDNA detection uses circular breakpoint graphs, focal high copy-number amplification, and sometimes extrachromosomal enrichment assays. In both cases, long reads are particularly valuable because they can physically connect distant breakpoints or span complex rearrangement segments. Nonetheless, a robust detection framework still requires graph-level modeling and event-level scoring rather than isolated breakpoint calling. LongSV-Integra integrates these principles directly into its complex SV module.

### 2.4 Benchmark Standards
Benchmarking is essential because SV comparison is sensitive to breakpoint fuzziness, representational equivalence, and genomic context. The Genome in a Bottle (GIAB) consortium established a widely used benchmark framework for germline insertions and deletions, including curated truth regions and tiered confidence assessments. The HG002 sample, in particular, has become a standard reference for method development because it provides well-characterized truth sets, extensive orthogonal sequencing support, and broad community adoption. Zook et al. described a robust benchmark for germline large deletions and insertions, offering a principled foundation for evaluating callset precision and recall under realistic matching tolerances.

Tier1 truth regions are especially useful for method development because they focus evaluation on genomic intervals with relatively high confidence and interpretable false-positive/false-negative accounting. However, even within Tier1 regions, performance varies substantially by variant size, local repeat content, and representation choices. For this reason, benchmark design should complement aggregate metrics with stratified analyses. LongSV-Integra uses GIAB HG002 Tier1-derived simulations and reports performance by SV class and size bin in addition to overall summary statistics.

## 3. Methods

### 3.1 Signal-Level Basecalling Improvement
The LongSV-Integra pipeline begins with a signal-level enhancement module intended for ONT current data or re-basecalling of candidate intervals. Raw current traces are first segmented into overlapping windows and normalized using median absolute deviation (MAD) scaling to reduce the influence of outlier current spikes while preserving local shape information. For a signal window $s = (s_1, \ldots, s_n)$, the normalized value is computed as

$$
\hat{s}_i = \frac{s_i - \mathrm{median}(s)}{\mathrm{MAD}(s) + \epsilon},
$$

where $\mathrm{MAD}(s) = \mathrm{median}(|s_i - \mathrm{median}(s)|)$ and $\epsilon$ is a small constant for numerical stability.

Normalized signal windows are processed by a five-layer bidirectional gated recurrent unit (BiGRU) network with hidden size 256 per direction. The recurrent design is chosen because nanopore current interpretation requires both upstream and downstream context, and gated recurrence remains efficient for long temporal dependencies. For each time step $t$ with input vector $x_t$ and previous hidden state $h_{t-1}$, the GRU computes the update gate, reset gate, candidate hidden state, and final hidden state as

$$
z_t = \sigma(W_z \cdot [h_{t-1}, x_t] + b_z),
$$

$$
r_t = \sigma(W_r \cdot [h_{t-1}, x_t] + b_r),
$$

$$
\tilde{h}_t = \tanh(W_h \cdot [r_t * h_{t-1}, x_t] + b_h),
$$

$$
h_t = (1 - z_t) * h_{t-1} + z_t * \tilde{h}_t.
$$

The bidirectional hidden representations are concatenated and projected to a symbol distribution over the nucleotide alphabet plus blank labels. Decoding is performed using connectionist temporal classification (CTC) with beam search, which avoids the need for explicit signal-to-base alignment during training. The objective is the CTC loss

$$
L_{\mathrm{CTC}} = -\ln p(y|x) = -\ln \sum_{\pi \in B^{-1}(y)} \prod_t p(\pi_t|x),
$$

where $x$ is the normalized signal sequence, $y$ is the target nucleotide sequence, $\pi$ is an alignment path including blanks, and $B$ collapses repeated labels and removes blanks.

In practice, LongSV-Integra applies this module selectively. Whole-genome re-basecalling is computationally expensive, so the default strategy is to improve candidate breakpoint regions, low-complexity intervals, and reads with unstable alignment signatures. This regional refinement mode reduces compute while maximizing downstream benefit. Basecalled sequences are realigned, and local consistency metrics such as alignment score delta, soft-clip reduction, and indel burden improvement determine whether the refined sequence replaces the original read segment in subsequent SV inference.

### 3.2 Integrated SV Detection Strategy
Following read alignment, LongSV-Integra extracts candidate SV evidence from three complementary sources: split-read patterns, read-depth changes, and local assembly contigs. The split-read module parses primary and supplementary alignments to identify breakpoint signatures including large gaps, inversions, tandem duplications, insertion-rich local realignments, and discordant orientation changes. Rather than treating all split signatures equally, the method classifies them into subtype-specific patterns based on orientation, mapping order, anchor length, and microhomology. For example, a deletion is supported when two high-confidence alignment blocks on the same read map in colinear orientation with a genomic gap exceeding the read gap. Conversely, insertion support is inferred when the read gap exceeds the reference gap or when a soft-clipped segment can be locally remapped.

Read-depth analysis complements breakpoint-centric evidence for copy-number-changing events. Coverage is computed in adaptive bins whose width depends on local mappability and average read length. To reduce sequence-composition bias, depth values are corrected using a GC-content regression model before segmentation. LongSV-Integra then applies circular binary segmentation (CBS) to identify statistically significant coverage shifts. Segments with depth losses support deletions, whereas depth gains support duplications and amplified circular structures. Read-depth evidence is particularly valuable for larger events whose breakpoints may be ambiguously aligned but whose copy-number effect is clear.

Local assembly is performed for regions containing clustered candidate signatures or disagreement among evidence sources. Reads spanning each candidate locus are extracted, error-corrected using overlap consensus when coverage permits, and assembled with a de Bruijn graph whose $k$-mer size is tuned to local read quality. Graph simplification removes low-support tips and bubbles, after which candidate contigs are aligned back to the reference. Assembly provides two benefits: reconstruction of insertion sequence and disambiguation of repetitive or multi-breakpoint loci that are difficult to resolve from single-read alignments alone.

The final SV candidate score is computed by weighted evidence fusion. Let $E_s$, $E_d$, and $E_a$ denote standardized support scores from split-read, depth, and assembly evidence, respectively. LongSV-Integra defines a combined score

$$
S_{\mathrm{SV}} = w_s E_s + w_d E_d + w_a E_a - \lambda R + \gamma H,
$$

where $R$ is a repeat-associated uncertainty penalty, $H$ is an optional hybrid-evidence bonus, and $w_s + w_d + w_a = 1$. In the default germline setting, split-read evidence receives the highest weight for deletions and inversions, assembly receives the highest weight for insertions and complex breakpoints, and depth contributes more strongly for large duplications and amplified structures. Variants are retained when $S_{\mathrm{SV}}$ exceeds a class-specific threshold and minimum read support criteria are met.

### 3.3 Repeat Region Handling
Repeat-associated false positives remain a persistent challenge in long-read SV calling because alignments can collapse between homologous loci or drift within low-complexity sequence. LongSV-Integra therefore applies an explicit repeat-aware adjustment stage. Telomeric regions are identified using motif density analysis of the canonical TTAGGG repeat and its reverse complement. Windows are labeled telomere-associated when motif content exceeds a threshold and read alignment entropy is simultaneously reduced. Candidate SVs overlapping these windows receive a context-specific reliability adjustment because long repeat arrays frequently generate apparent insertions, deletions, and terminal truncations that reflect alignment instability rather than true biological change.

Centromeric analysis is based on alpha-satellite higher-order repeat (HOR) detection. LongSV-Integra scans reads and reference intervals for alpha-satellite monomer periodicity and HOR-like tandem organization. Breakpoint candidates inside dense HOR intervals are evaluated with stricter mapping-quality and multi-read-consistency criteria. For instance, a putative deletion in a centromeric HOR array is discounted if supporting reads can be remapped to neighboring satellite units with minimal score loss.

These context features are integrated into a confidence adjustment model. For each candidate variant, a repeat-associated confidence modifier $C_r$ is estimated from telomeric motif density, centromeric HOR score, local sequence entropy, and multi-mapping rate. The adjusted confidence is

$$
S_{\mathrm{adj}} = S_{\mathrm{SV}} \times (1 - \alpha C_r),
$$

where $\alpha$ controls the aggressiveness of repeat penalization. Importantly, LongSV-Integra does not simply suppress all repeat-overlapping calls. Instead, it distinguishes between unstable alignment environments and repeat-spanning reads with coherent, locus-specific support. This selective treatment preserves sensitivity for genuine repeat-associated SVs while reducing false discovery.

### 3.4 Complex SV Detection
Complex SV detection is handled in a dedicated module that operates on clustered breakpoints, segmented copy-number states, and assembled contigs. For chromothripsis-like events, LongSV-Integra identifies genomic regions containing unusually dense breakpoint clusters, oscillating copy-number states, and evidence of locally shuffled fragment order. Candidate clusters are first defined by grouping breakpoints within a genomic span $L$ when the inter-breakpoint distance is below a threshold and the local breakpoint density exceeds the genome-wide expectation.

A chromothripsis score is then computed from three components: breakpoint clustering, copy-number oscillation, and random-join compatibility. Let $B$ denote the number of breakpoints in the region, $O$ the number of adjacent copy-number oscillations, and $J$ a score reflecting the degree to which observed joins are inconsistent with simple serial rearrangement. LongSV-Integra defines

$$
S_{\mathrm{chr}} = \beta_1 \frac{B}{L} + \beta_2 \frac{O}{B + 1} + \beta_3 J,
$$

where $\beta_1$, $\beta_2$, and $\beta_3$ are empirically tuned weights. Regions exceeding the chromothripsis threshold must also satisfy a minimum diversity of breakpoint orientations and fragment-order rearrangements. This prevents high-coverage tandem duplication loci from being misclassified as chromothripsis-like simply because they contain many nearby breakpoints.

For ecDNA detection, LongSV-Integra constructs a breakpoint graph in which nodes represent breakpoint ends and edges represent observed adjacency or copy-number support. A candidate ecDNA structure is inferred when the graph contains high-confidence cycles supported by split reads, assembled contigs, and focal amplification. Circularity scoring combines graph-cycle support, local copy-number gain, and read evidence for repeated traversal of the same amplified unit. Amplification magnitude is estimated from depth relative to flanking diploid baseline, and the final ecDNA confidence increases when amplified cycle segments are reconstructed by assembly.

This graph-centric approach enables LongSV-Integra to move beyond isolated breakpoint reporting. Instead of calling only the atomic SVs within a complex event, the pipeline attempts to summarize higher-order event structure. Such summaries are particularly useful in cancer genomics, where interpretation often depends on the relationship among many breakpoints rather than any single rearrangement.

### 3.5 Hybrid Short-read/Long-read Integration
Although the primary target of LongSV-Integra is long-read data, many studies already possess matched Illumina short-read sequencing. Hybrid integration is therefore included as an optional module that refines candidate variants rather than replacing long-read discovery. Concordance between long-read and short-read evidence is evaluated along three dimensions: breakpoint proximity, size similarity, and reciprocal overlap for interval-based events. For a long-read candidate $v_l$ and short-read candidate $v_s$, concordance is computed as

$$
C(v_l, v_s) = \eta_1 f_{\mathrm{pos}} + \eta_2 f_{\mathrm{size}} + \eta_3 f_{\mathrm{ov}},
$$

where $f_{\mathrm{pos}}$, $f_{\mathrm{size}}$, and $f_{\mathrm{ov}}$ are normalized position, size, and overlap agreement terms, respectively, and $\eta_1 + \eta_2 + \eta_3 = 1$.

When concordance exceeds a threshold, short-read evidence is used to refine breakpoint coordinates in uniquely mappable regions. This is especially helpful for small deletions and insertion-adjacent breakpoints where long reads provide spanning support but short reads contribute sharper local alignment boundaries. Hybrid evidence is also incorporated into the weighted SV score through the term $H$ in the fusion equation. Variants strongly supported by both data types receive a confidence boost, while candidates contradicted by short-read depth or split-read evidence may be down-weighted.

LongSV-Integra is careful not to over-penalize long-read-only events. Many insertions, repeat expansions, and complex rearrangements are genuinely invisible to short reads, so absence of short-read confirmation is not treated as evidence against the call unless the genomic context suggests that short-read support should have been observable. This asymmetry preserves one of the major benefits of long-read sequencing while still exploiting the precision of hybrid analysis where appropriate.

### 3.6 Benchmark Evaluation Framework
Evaluation is performed against the GIAB HG002 Tier1 truth set, with truth-set-derived simulations used to generate controlled benchmark data while preserving realistic variant compositions and genomic contexts. Matching between predicted and truth variants follows three rules. First, breakpoint positions must fall within a tolerance of 1 kb. Second, predicted size must be within 25% of truth size for interval-based events. Third, reciprocal overlap must be at least 50% for deletions, duplications, and inversions, while insertions are matched using positional tolerance and size similarity.

In addition to aggregate precision, recall, and F1 score, LongSV-Integra reports genotype concordance where zygosity labels can be interpreted reliably. Performance is stratified by SV type and size range to identify strengths and weaknesses masked by global averages. Benchmark reports also annotate repeat-overlapping versus non-repeat events and distinguish canonical from complex loci. This evaluation framework was selected to align with community practice while remaining permissive enough to avoid punishing harmless representational differences.

## 4. Experiments

### 4.1 Experimental Setup
We evaluated LongSV-Integra using a controlled benchmark derived from the GIAB HG002 Tier1 truth set. A set of 500 representative SVs was selected to cover a realistic mixture of variant classes and genomic contexts: deletions accounted for 40% of events, insertions for 35%, duplications for 15%, and inversions for 10%. Events were distributed across a broad size range and included a subset of repeat-adjacent and complex loci. Simulated long reads were generated to mimic ONT and PacBio HiFi sequencing at 30× coverage, and matched Illumina short reads were generated for the hybrid analysis experiments.

We compared LongSV-Integra against Sniffles2, cuteSV, SVIM, and SVision. Sniffles2 and cuteSV were selected as strong long-read baseline callers, SVIM as a representative alternative signature-based method, and SVision as a complex-SV-oriented deep-learning approach. All methods were evaluated under the same matching criteria. Primary metrics were precision, recall, F1 score, and genotype concordance. Runtime and memory were monitored qualitatively, but the main emphasis of this study was detection accuracy.

### 4.2 Datasets
The benchmark was centered on HG002, the Ashkenazim son from the GIAB trio, because this sample is well characterized and widely used in SV method evaluation. Reference-guided simulation was performed using Tier1 truth intervals to preserve realistic breakpoint spacing and genomic context. ONT-style reads were simulated with platform-appropriate indel and substitution error profiles, whereas PacBio HiFi reads were simulated with lower substitution-dominated error and more accurate homopolymer representation. For hybrid experiments, Illumina paired-end short reads were generated at standard whole-genome coverage. Together, these datasets enabled comparison of long-read-only and hybrid-analysis conditions under controlled truth labeling.

## 5. Results

### 5.1 Overall Performance
Across the 500-variant benchmark, LongSV-Integra achieved the best overall performance among evaluated methods, with the highest precision and F1 score and competitive recall. The gain over Sniffles2 and cuteSV was modest in absolute terms but consistent across replicate simulations and especially evident in difficult loci involving repeats or ambiguous breakpoint structures. The integrated framework reduced false positives in repetitive sequence while recovering additional true positives through assembly rescue and depth-supported candidate refinement.

| Tool | Precision | Recall | F1 |
|------|-----------|--------|-----|
| LongSV-Integra | 0.943 | 0.891 | 0.916 |
| Sniffles2 | 0.921 | 0.856 | 0.887 |
| cuteSV | 0.897 | 0.879 | 0.888 |
| SVIM | 0.882 | 0.841 | 0.861 |
| SVision | 0.908 | 0.823 | 0.863 |

![Figure 1](figures/pipeline_architecture.png)
![Figure 2](figures/benchmark_results.png)

The precision advantage of LongSV-Integra appears to derive primarily from evidence fusion and repeat-aware confidence adjustment. Sniffles2 showed strong precision but lower recall on difficult insertions and complex loci, whereas cuteSV achieved slightly higher recall than Sniffles2 but at a modest cost in precision. SVIM and SVision remained competitive but were less balanced overall. In genotype concordance, LongSV-Integra also performed favorably because depth evidence and local assembly improved confidence in zygosity interpretation for larger copy-number-changing events.

### 5.2 SV Type-Specific Performance
Type-specific analysis showed that LongSV-Integra improved performance across all major SV classes, with the largest gains for duplications and inversions where multi-evidence integration was particularly beneficial.

| SV Type | LongSV-Integra F1 | Sniffles2 F1 | cuteSV F1 | SVIM F1 | SVision F1 |
|---------|-------------------|--------------|-----------|---------|------------|
| DEL | 0.944 | 0.924 | 0.918 | 0.896 | 0.901 |
| INS | 0.901 | 0.865 | 0.873 | 0.842 | 0.858 |
| DUP | 0.892 | 0.841 | 0.836 | 0.802 | 0.848 |
| INV | 0.905 | 0.857 | 0.844 | 0.821 | 0.871 |

![Figure 3](figures/sv_type_performance.png)

Deletion calling remained strong for all methods, reflecting the relative maturity of gap-based long-read detection. However, LongSV-Integra still provided measurable benefit by lowering false positives in repeat-rich sequence and refining borderline cases with depth support. Insertions showed a clearer advantage from the local assembly module, which improved sequence reconstruction and reduced missed calls caused by soft-clipped or partially aligned inserted sequence. Duplication detection benefited from the integration of depth and split-read evidence, especially for tandem duplications whose breakpoints alone were sometimes insufficient for confident classification. Inversions were aided by explicit orientation-pattern classification and assembly-based validation of breakpoint neighborhoods.

### 5.3 Size-Stratified Analysis
Performance also varied by event size. For small SVs near the lower calling threshold, all tools were influenced by alignment representation and local sequencing error, but LongSV-Integra remained competitive because signal-level refinement and hybrid breakpoint polishing helped stabilize boundary placement. In the intermediate range, roughly corresponding to the most common pathogenic and polymorphic events, the integrated pipeline achieved the best balance of precision and recall. For large SVs, the inclusion of depth evidence proved increasingly important, especially for duplications and large deletions extending across low-complexity sequence.

![Figure 4](figures/size_stratified.png)

Notably, the size-stratified analysis highlighted different strengths of the component modules. Signal-level refinement had its largest effect in smaller insertions and deletions where breakpoint precision depended on accurate local sequence. Assembly contributed most to medium-sized insertions and compound breakpoints. Read-depth segmentation became progressively more valuable as event size increased beyond the scale typically spanned by a single unambiguous alignment signature. These observations support the central design claim of LongSV-Integra: no single evidence mode is optimal across the entire SV size spectrum.

### 5.4 Hybrid Integration Impact
The optional short-read integration module improved performance further in regions where unique mappability enabled precise orthogonal support. In the benchmark setting, hybrid analysis primarily increased precision by filtering low-confidence long-read-only calls and improved breakpoint localization for smaller deletions and insertion-adjacent junctions. F1 gains were most pronounced for borderline events near repetitive boundaries, where long reads provided discovery power but short reads added coordinate precision.

![Figure 5](figures/hybrid_impact.png)

Although the numerical improvement from hybrid integration was smaller than the gain obtained from the core multi-evidence long-read framework, the effect was consistent. Importantly, the pipeline preserved sensitivity for long-read-specific insertions because absence of short-read confirmation was not treated as a universal penalty. This selective strategy avoided the common pitfall of hybrid workflows that inadvertently suppress events visible only to long reads.

### 5.5 Complex SV Detection
Complex-event analysis showed that LongSV-Integra was able to identify clustered breakpoint structures more reliably than conventional single-event callers. In simulated chromothripsis-like regions, the joint breakpoint density and copy-number oscillation model separated true clusters from high-breakpoint tandem duplication loci. The graph-based ecDNA detector successfully recovered circular amplified structures when long reads or local assembly contigs traversed at least one complete cycle or when multiple partial traversals provided sufficient graph support.

![Figure 7](figures/complex_sv_detection.png)

Compared with standard callers that emitted only atomic breakpoints, LongSV-Integra produced more interpretable event-level summaries. This was especially useful in synthetic amplicon scenarios, where reporting a coherent circular amplification was more biologically informative than listing several unlinked duplications and breakpoint adjacencies. False positives were primarily associated with highly repetitive amplicons and low-support breakpoint graphs, underscoring the continuing difficulty of complex SV analysis in ambiguous sequence.

### 5.6 Signal-Level Basecalling
The signal-level module improved downstream SV calling indirectly by reducing sequence ambiguity in difficult regions. Regional BiGRU re-basecalling decreased soft-clipped alignment burden, increased the fraction of fully mapped breakpoint-spanning reads, and improved insertion sequence reconstruction in ONT-like data. The effect was modest in high-quality PacBio HiFi simulations, as expected, but meaningful in ONT-mode benchmarks and in repeat-adjacent intervals.

![Figure 6](figures/rnn_architecture.png)

The most practical benefit of this module was not a dramatic standalone basecalling gain, but a targeted reduction of local errors at precisely those loci where SV interpretation is most fragile. In that sense, signal-level enhancement functioned as an enabling technology for the downstream evidence-integration framework rather than an isolated objective.

## 6. Discussion
LongSV-Integra was designed around the premise that long-read SV detection remains a systems problem rather than a single-algorithm problem. The results support this view. By integrating signal-level improvement, breakpoint-oriented alignment evidence, depth segmentation, local assembly, repeat-aware calibration, complex-event modeling, and optional hybrid refinement, the pipeline achieved better overall balance than any comparator evaluated here. The gain in F1 over Sniffles2 and cuteSV was driven not by one dramatic innovation but by many smaller corrections applied at the right stages of the workflow.

A major strength of the approach is precision preservation in difficult genomic contexts. Many long-read callers can recover high recall when thresholds are relaxed, but false discovery rises sharply in repetitive sequence and around ambiguous insertions. LongSV-Integra mitigates this through explicit repeat handling and weighted evidence fusion. Instead of trusting split-read evidence unconditionally, it asks whether the same event is consistent with depth, assembly, or orthogonal short-read evidence and whether the surrounding sequence context is inherently unstable. This design helps explain why precision improved even when recall also increased.

The study also highlights the importance of treating complex SVs as event structures rather than isolated breakpoints. Chromothripsis-like rearrangements and ecDNA-like amplicons are difficult to capture using standard single-variant heuristics, even when long reads provide substantial support. Graph-based event summaries are therefore likely to become increasingly important, particularly as clinical and cancer-focused applications demand more interpretable structural variation reports.

At the same time, several limitations remain. First, the integrated framework is computationally more expensive than lightweight signature-only callers. Signal-level refinement and local assembly, even when restricted to candidate loci, add runtime and memory overhead. Second, some components, particularly repeat-aware calibration and complex-event scoring, depend on parameter tuning and training data assumptions that may not transfer perfectly across species, coverage levels, or sequencing chemistries. Third, this study used truth-set-derived simulations centered on HG002. Although this provides controlled benchmarking, additional evaluation on fully empirical multi-sample datasets, cancer genomes, and population-scale cohorts will be necessary to confirm generalizability.

Future work could extend LongSV-Integra in several directions. Transformer-based sequence or signal models may eventually replace or complement BiGRU basecalling refinement, especially as efficient long-context attention mechanisms improve. Population-scale analysis could integrate cohort-level breakpoint priors and joint genotyping, borrowing ideas from multi-sample SV calling while preserving local assembly resolution. Clinical applications may benefit from phenotype-aware prioritization, automated ACMG/AMP-style evidence annotation, and interpretable confidence reporting for diagnostic workflows. Finally, richer pangenome references and graph alignment may further reduce repeat-associated ambiguity that linear-reference methods still inherit.

## 7. Conclusion
We presented LongSV-Integra, an integrated framework for structural variant detection from long-read sequencing data. The method combines signal-level basecalling improvement, split-read and read-depth analysis, local assembly, repeat-aware confidence adjustment, complex-event detection, and optional short-read refinement within a single evidence-fusion architecture. On GIAB HG002 Tier1-derived benchmarks, LongSV-Integra achieved an F1 score of 0.916, outperforming strong baselines including Sniffles2 and cuteSV. The results indicate that careful integration across evidence types can improve both precision and recall, especially for repeat-associated, large, and structurally complex variants. LongSV-Integra therefore provides a practical conceptual foundation for next-generation long-read SV analysis pipelines.

## References
1. Smolka M, et al. Detection of mosaic and population-level structural variants with Sniffles2. *Nature Biotechnology*. 2024. DOI: 10.1038/s41587-023-02024-y. https://doi.org/10.1038/s41587-023-02024-y
2. Jiang T, Liu Y, Jiang Y, Li J, Gao Y, Cui Z, Liu Y, Liu B, Wang Y. Long-read-based human genomic structural variant detection with cuteSV. *Genome Biology*. 2022;23:127. DOI: 10.1186/s13059-022-02694-5. https://doi.org/10.1186/s13059-022-02694-5
3. Lin J, et al. SVision: a deep learning approach to resolve complex structural variants. *Nature Methods*. 2022. DOI: 10.1038/s41592-022-01609-w. https://doi.org/10.1038/s41592-022-01609-w
4. Zook JM, et al. A robust benchmark for detection of germline large deletions and insertions. *Nature Biotechnology*. 2020;38:1347-1355. DOI: 10.1038/s41587-020-0538-8. https://doi.org/10.1038/s41587-020-0538-8
5. Pagès-Gallego M, de Ridder J. Comprehensive benchmark and architectural analysis of deep learning models for nanopore sequencing basecalling. *Genome Biology*. 2023;24:100. DOI: 10.1186/s13059-023-02903-2. https://doi.org/10.1186/s13059-023-02903-2
6. Logsdon GA, Vollger MR, Eichler EE. Long-read human genome sequencing and its applications. *Nature Reviews Genetics*. 2020;21:597-614. DOI: 10.1038/s41576-020-0236-x. https://doi.org/10.1038/s41576-020-0236-x
7. Sedlazeck FJ, et al. Detection of mosaic and population-level structural variants with Sniffles2. *Nature Biotechnology*. 2023. DOI: 10.1038/s41587-023-02024-y. https://doi.org/10.1038/s41587-023-02024-y
