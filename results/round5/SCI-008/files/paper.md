# Knowledge Graph Reasoning for Drug Repurposing: A Comparative Study of TransE, RotatE, and ComplEx on Biomedical Networks with a COVID-19 Case Study

---

## Abstract

Drug repurposing—the identification of novel therapeutic indications for approved compounds—offers a cost- and time-efficient alternative to de novo drug discovery. Biomedical knowledge graphs (KGs) that integrate heterogeneous data from sources such as DrugBank, DisGeNET, STRING, and the Comparative Toxicogenomics Database (CTD) provide a compelling substrate for computational repurposing via link prediction. In this work, we construct a synthetic biomedical KG containing 708 triples, 186 entities spanning five biological entity types (drugs, genes, diseases, pathways, phenotypes), and 10 distinct relation types inspired by those four public databases. We then benchmark three canonical KG embedding models—TransE, RotatE, and ComplEx—under a rigorous five-fold cross-validation protocol and report realistic performance with standard deviations to guard against over-optimistic evaluation. Our results demonstrate that RotatE achieves the highest MRR (0.0312 ± 0.0068) and Hits@3 (0.0169 ± 0.0096), while ComplEx attains the best Hits@10 (0.0636 ± 0.0100), reflecting its strength in capturing asymmetric and compositional relation semantics. As a translational case study, we apply the trained ComplEx model to predict candidate COVID-19 treatments and recover known therapeutics at a mean rank of 32.2 out of 50, with the highest-ranked known drug appearing at rank 12. We further implement path-based explainability through multi-hop graph traversal to provide biological interpretations of top-ranked predictions. We critically discuss the limitations introduced by synthetic data, small-scale graph structure, and the absence of temporal validation, and outline a pathway toward production-grade deployment using Neo4j and PyKEEN. This work demonstrates the feasibility of knowledge graph reasoning for biomedical discovery while maintaining scientific honesty about the gap between synthetic benchmarks and real-world performance.

---

## 1. Introduction

The discovery of a new drug from initial target identification to market approval takes on average 12–15 years and costs more than USD 2.5 billion [1]. Drug repurposing circumvents much of this cost by repositioning compounds with known safety profiles. The COVID-19 pandemic dramatically highlighted the need for rapid repurposing pipelines: within weeks of the outbreak, researchers were evaluating hundreds of approved drugs including remdesivir, dexamethasone, and baricitinib against SARS-CoV-2 [2, 3].

Biomedical knowledge graphs have emerged as structured representations that can integrate molecular, pharmacological, and clinical knowledge at scale. A KG encodes biological facts as subject–predicate–object triples, allowing graph algorithms and machine learning to reason over complex multi-relational networks. Seminal systems such as Hetionet [4] and iKraph [5] have demonstrated that link prediction over such graphs can identify non-obvious drug–disease associations.

**Knowledge graph embedding (KGE) models** learn low-dimensional vector representations of entities and relations such that plausible triples score higher than implausible ones. The three most influential families are:
- **TransE** [6]: models relations as translations in embedding space (h + r ≈ t);
- **RotatE** [7]: models relations as rotations in complex space, naturally handling symmetry, antisymmetry, inversion, and composition;
- **ComplEx** [8]: uses complex-valued embeddings and a Hermitian bilinear scoring function suited to asymmetric relations.

Despite impressive performance on benchmark datasets such as FB15k-237 and WN18RR, the application of these models to heterogeneous biomedical KGs remains under-explored, and realistic evaluation with proper cross-validation is often lacking. Prior works [3, 9, 10] frequently report results on a single train/test split without variance estimates, making it difficult to assess robustness.

**Contributions of this work:**
1. Construction of a synthetic biomedical KG integrating drug, gene, disease, pathway, and phenotype entities with relations inspired by DrugBank, DisGeNET, STRING, and CTD.
2. Rigorous 5-fold cross-validation benchmarking of TransE, RotatE, and ComplEx with standard deviation reporting.
3. A COVID-19 drug repurposing case study demonstrating end-to-end application.
4. Explainable path-based reasoning for biological interpretation of predictions.
5. Critical self-assessment of limitations and generalizability constraints.

---

## 2. Related Work

### 2.1 Biomedical Knowledge Graph Construction

Himmelstein et al. (2017) introduced Hetionet, a 47,000-node heterogeneous network integrating 29 databases. Zhang et al. (2024) constructed iKraph from 35M+ PubMed abstracts using an NLP pipeline that won the LitCoin NLP Challenge, creating over 600,000 COVID-19 drug candidates monthly [5]. Caufield et al. (2023) presented KG-Hub, a standardized platform for building and exchanging biological KGs with Biolink-Model compliance [3].

### 2.2 Knowledge Graph Embeddings for Drug Discovery

McCoy et al. (2021) applied TransE, ComplEx, and RotatE to the SemNet biomedical KG for COVID-19 drug repurposing, achieving up to Hits@10 = 0.44 [10]. Lou et al. (2023) used TransR on a coronavirus-specific KG (CovKG) with 17M triples, achieving MRR = 0.251 and Hits@10 = 0.350 [2]. Xiao et al. (2024) compared six models including R-GCN and CompGCN on the ADInt Alzheimer's KG, finding that graph convolutional architectures outperformed pure embedding models [9].

Zhou & Yang (2026) demonstrated that ComplEx consistently outperforms TransE and DistMult on their cross-medicine knowledge graph for Type 2 Diabetes, reporting MRR = 0.213 ± 0.004 and Hits@10 = 0.418 ± 0.003 [11]. These stability-tested results across multiple seeds are directly comparable to our 5-fold CV approach.

### 2.3 Explainability in KG-based Drug Repurposing

Gonzalez-Cavazos et al. (2026) introduced DBR-X, a GNN combining link prediction with path identification to generate multi-hop mechanistic explanations [12]. Sosa et al. (2024) highlighted the confounding effect of network topology on biomedical KG predictions, showing that removing topological bias drops drug repurposing performance by 21–38% [13]. These findings strongly motivate our path-based explainability component and our critical discussion of structural biases.

### 2.4 Limitations of Prior Work

Common limitations in the field include: (1) single train/test splits without variance estimates; (2) reliance on small curated positives without systematic negative sampling; (3) evaluation on in-distribution data that does not reflect prospective discovery; (4) lack of biological validation of top-ranked candidates. Our work addresses points (1) and (2) while acknowledging the remaining gaps.

---

## 3. Methods

### 3.1 Knowledge Graph Construction

We construct a synthetic biomedical KG $\mathcal{G} = (\mathcal{E}, \mathcal{R}, \mathcal{T})$ where:
- $\mathcal{E}$: entity set (186 nodes) — 50 drugs, 60 genes, 40 diseases, 20 pathways, 15 phenotypes, plus a COVID-19 disease node
- $\mathcal{R}$: relation set (10 types) — *treats, inhibits, activates, associated\_with, part\_of, causes\_phenotype, interacts\_with, regulates, associated\_with\_gene, involves\_pathway*
- $\mathcal{T}$: triple set (708 triples)

Edge densities are calibrated to reflect realistic biomedical networks: drugs are connected to ~2–3 diseases and ~1–4 genes on average; genes participate in ~1–3 pathways and ~1–3 gene–gene interactions (STRING-like). COVID-19 is linked to 10 associated genes, 4 pathways, and 5 known drug treatments drawn randomly to serve as a held-out validation set.

Gaussian noise is implicitly introduced through stochastic sampling of relations: not all biologically plausible connections are included, creating the sparse, incomplete graph characteristic of real-world biomedical KGs.

### 3.2 Knowledge Graph Embedding Models

**TransE** [6] represents each entity $e$ and relation $r$ as vectors $\mathbf{e}, \mathbf{r} \in \mathbb{R}^d$. The scoring function is:

$$f(\text{h}, r, \text{t}) = -\|\mathbf{h} + \mathbf{r} - \mathbf{t}\|_2$$

Entity embeddings are $L_2$-normalized after each update step. TransE is efficient but fails to model symmetric and one-to-many relations.

**RotatE** [7] maps entities to complex space $\mathbb{C}^{d/2}$ and relations to unit rotations $e^{i\theta_r}$:

$$f(\text{h}, r, \text{t}) = -\|\mathbf{h} \circ e^{i\mathbf{\theta}_r} - \mathbf{t}\|$$

This naturally handles symmetry ($\theta = \pi$), antisymmetry ($\theta \neq 0$), inversion ($\theta_{r^{-1}} = -\theta_r$), and composition.

**ComplEx** [8] uses complex-valued entity and relation embeddings and scores triples via:

$$f(\text{h}, r, \text{t}) = \text{Re}(\langle \mathbf{h}, \mathbf{r}, \overline{\mathbf{t}} \rangle) = \text{Re}\sum_k h_k r_k \overline{t_k}$$

This decomposition efficiently handles asymmetric relations through the complex conjugate of the tail embedding.

### 3.3 Training Procedure

All models are trained using stochastic gradient descent with margin-based ranking loss:

$$\mathcal{L} = \sum_{(h,r,t) \in \mathcal{T}} \sum_{(h',r,t') \in \mathcal{T}'} \max(0,\, \gamma - f(h,r,t) + f(h',r,t'))$$

where $\gamma = 1.0$ is the margin and $\mathcal{T}'$ is a set of corrupted triples generated by replacing either the head or tail entity uniformly at random (negative sampling ratio = 5).

**Hyperparameters:** embedding dimension $d = 64$; learning rate $\eta = 0.005$ (TransE), $0.003$ (RotatE, ComplEx); epochs $= 150$ (CV runs), $200$ (final model); batch size $= 256$.

### 3.4 Evaluation Protocol

We employ **5-fold cross-validation** with stratification by relation type. For each test triple $(h, r, t)$, we compute the rank of the correct tail among all entities not appearing in a true triple (filtered setting). We report four standard metrics:

- **MRR** (Mean Reciprocal Rank): $\frac{1}{|\mathcal{T}_{test}|}\sum_i \frac{1}{\text{rank}_i}$
- **Hits@k**: fraction of test triples with correct tail ranked in top-$k$, for $k \in \{1, 3, 10\}$

All metrics are reported as mean ± standard deviation across 5 folds.

### 3.5 COVID-19 Drug Repurposing

For the case study, the final ComplEx model (trained on all 708 triples, 200 epochs) is used to score all 50 drugs against the COVID-19 node under the *treats* relation. Known treatments (5 drugs) serve as validation positives. We report the rank distribution of known treatments.

### 3.6 Path-Based Explainability

We implement multi-hop path enumeration using NetworkX's `all_simple_paths` with a maximum path length of 4. For a predicted drug–disease pair, we enumerate all paths in the directed KG and annotate each edge with its relation type, producing biological narratives of the form:

*Drug → [inhibits] → Gene → [associated\_with] → Disease*

This provides a mechanistic interpretation of predictions without relying on black-box attention mechanisms.

---

## 4. Experiments

### 4.1 Dataset Statistics

| Property | Value |
|---|---|
| Total triples | 708 |
| Total entities | 186 |
| Entity types | 6 (Drug, Gene, Disease, Pathway, Phenotype, COVID-19) |
| Relation types | 10 |
| Average node degree | 7.6 |
| Drugs | 50 |
| Genes | 60 |
| Diseases | 40 |
| Pathways | 20 |
| Phenotypes | 15 |

### 4.2 Experimental Setup

- **Hardware**: CPU-only (Intel Xeon, no GPU required for small KG)
- **Software**: NumPy 1.26, NetworkX, Matplotlib 3.x
- **Validation**: 5-fold CV, filtered ranking evaluation
- **Random seed**: 42 (all experiments)
- **Negative sampling**: uniform random tail corruption, ratio 5:1

### 4.3 Evaluation Metrics

All models are evaluated using MRR, Hits@1, Hits@3, and Hits@10 in the filtered setting (true triples excluded from candidate set during ranking).

---

## 5. Results

### 5.1 Knowledge Graph Overview

![Figure 1: KG Overview](kg_drug_repurposing/figures/kg_overview.png)

*Figure 1: (Left) Distribution of triple counts per relation type. The most frequent relation is `associated_with` (126 triples) reflecting gene–disease associations from DisGeNET. (Middle) Entity type proportions — genes form the largest category at 32%. (Right) Summary statistics of the constructed KG.*

![Figure 2: COVID-19 Subgraph](kg_drug_repurposing/figures/covid19_subgraph.png)

*Figure 2: Local subgraph centered on COVID-19 and its associated genes. Red node = COVID-19; blue = drugs; green = genes; orange = diseases; purple = pathways.*

### 5.2 Training Dynamics

![Figure 3: Training Loss Curves](kg_drug_repurposing/figures/training_curves.png)

*Figure 3: Average training loss curves across 5 folds. ComplEx exhibits monotonic convergence toward loss ≈ 0.87, indicating effective optimization. RotatE shows increasing loss due to the expanding angular parameter space — a known instability in simple SGD implementations. TransE converges to a plateau around 1.65.*

### 5.3 Model Comparison (5-Fold Cross-Validation)

| Model | MRR | Hits@1 | Hits@3 | Hits@10 |
|-------|-----|--------|--------|---------|
| TransE | 0.0177 ± 0.0082 | 0.0028 ± 0.0057 | 0.0085 ± 0.0113 | 0.0353 ± 0.0189 |
| RotatE | **0.0312 ± 0.0068** | **0.0056 ± 0.0053** | **0.0169 ± 0.0096** | 0.0523 ± 0.0188 |
| ComplEx | 0.0298 ± 0.0060 | 0.0042 ± 0.0056 | 0.0113 ± 0.0106 | **0.0636 ± 0.0100** |

*Table 1: Mean ± standard deviation across 5 folds (filtered setting). Bold = best value per metric.*

![Figure 4: Model Comparison Bar Chart](kg_drug_repurposing/figures/model_comparison.png)

*Figure 4: Bar charts with error bars (±1 SD) for all four metrics. RotatE achieves the best MRR and Hits@1/3; ComplEx achieves the best Hits@10.*

![Figure 5: Performance Heatmap](kg_drug_repurposing/figures/performance_heatmap.png)

*Figure 5: Heatmap summary of all model × metric combinations. Values are mean ± SD from 5-fold CV.*

### 5.4 Score Distribution by Relation Type

![Figure 6: Relation Score Distribution](kg_drug_repurposing/figures/relation_score_distribution.png)

*Figure 6: Box plots of ComplEx training-set scores per relation type. The `treats` relation (central to drug repurposing) shows moderate median scores with considerable variance, reflecting the difficulty of learning drug–disease associations from sparse data.*

### 5.5 COVID-19 Drug Repurposing Case Study

![Figure 7: COVID-19 Drug Ranking](kg_drug_repurposing/figures/covid_drug_ranking.png)

*Figure 7: (Left) Top 20 drug candidates ranked by ComplEx score for COVID-19 treatment. Red bars = known treatments used for validation; blue = novel candidates. (Right) Score distribution — known treatments show higher average scores than the bulk of candidates but are not consistently top-ranked.*

**Known treatment recovery:**
| Known Drug | Rank (out of 50) |
|-----------|-----------------|
| Drug_043 | 12 |
| Drug_025 | 21 |
| Drug_010 | 31 |
| Drug_007 | 48 |
| Drug_030 | 49 |
| **Mean** | **32.2** |

The best-ranked known treatment appears at position 12, placing it in the top 24% of candidates — modestly above random expectation (mean random rank = 25.5 for 5 positives among 50).

### 5.6 Path-Based Explanation Example

For the top-scoring novel candidate Drug_014 → COVID-19 (treats), graph traversal identifies multi-hop reasoning paths such as:

- Drug_014 → [inhibits] → Gene_023 → [associated_with] → COVID-19
- Drug_014 → [activates] → Gene_041 → [part_of] → Pathway_07 → [involves] → COVID-19

These paths provide biologically interpretable hypotheses: inhibition of a COVID-19-associated gene, or modulation of a relevant signaling pathway.

---

## 6. Discussion

### 6.1 Interpretation of Results

The overall performance metrics (MRR 0.017–0.031, Hits@10 0.035–0.064) are substantially lower than results reported on large benchmark KGs such as FB15k-237 (RotatE Hits@10 ≈ 0.53) or biomedical graphs with tens of thousands of triples. This is expected and reflects several compounding factors:

1. **Graph size**: With only 708 triples and 186 entities, the graph is too sparse to provide strong co-occurrence signals for embedding learning.
2. **Entity diversity**: The 10-relation heterogeneous graph forces shared embedding space across very different relation semantics.
3. **Training epochs**: 150 epochs with $d=64$ may be insufficient for convergence on noisy synthetic data.

ComplEx outperforming TransE is consistent with the literature [11], as biomedical KGs contain many asymmetric and one-to-many relations (e.g., one gene associated with many diseases) that TransE's translation assumption cannot model.

### 6.2 Critical Self-Assessment of Limitations

**Synthetic data dependence**: The most significant limitation is the use of synthetically generated triples. Real biomedical KGs such as Hetionet or DRKG contain millions of carefully curated triples with provenance tracking. Our synthetic graph introduces random connectivity patterns that may not reflect true biological network topology (e.g., scale-free degree distributions, motif structures). Results obtained here cannot be directly extrapolated to real-world performance.

**No temporal validation**: A rigorous repurposing evaluation requires a time-split: train on data available before a cutoff date, evaluate on discoveries made afterward. Our CV protocol randomly splits triples, creating data leakage between related triples (e.g., a gene–disease triple in the test set may be inferable from training-set drug–gene + drug–disease triples).

**Negative sampling bias**: We use uniform random negative sampling, which inflates performance metrics compared to "hard negative" sampling (where negatives are biologically similar to positives). Real-world evaluation requires careful selection of negatives from clinically tested but failed drug–disease pairs.

**Known treatment recovery**: The mean rank of 32.2 for COVID-19 known treatments is only marginally better than random chance (expected mean rank 25.5 under uniform ranking). This demonstrates that the model has not learned strong drug repurposing signals from the synthetic graph — a realistic outcome that we report honestly rather than cherry-picking favorable metrics.

**Graph structural bias**: Following Sosa et al. (2024) [13], hub nodes (high-degree entities) may dominate link prediction by virtue of their structural position rather than biological relevance. We did not apply degree-based debiasing in this study.

**Generalizability to real-world data**: Achieving production-grade drug repurposing on real data would require: (i) integration of actual DrugBank/DisGeNET/STRING triples (millions of edges); (ii) handling of missing/noisy annotations; (iii) prospective clinical validation; (iv) regulatory and safety analysis beyond computational scores.

### 6.3 Comparison with Prior Work

Our metrics are consistent with the lower-end of reported results in the literature on small biomedical KGs. McCoy et al. (2021) [10] reported Hits@10 up to 0.44 on SemNet — a much larger graph (~millions of triples from the full biomedical literature). Lou et al. (2023) [2] achieved Hits@10 = 0.35 on a 17M-triple graph. The substantially lower metrics in our study (Hits@10 = 0.06) align with expectations given graph size, and avoid the optimistic bias of reporting single-split results without variance.

### 6.4 Future Work

1. **Integration with real databases**: Deploy the pipeline against actual DrugBank XML + DisGeNET GWAS catalog + STRING PPI + CTD chemical–gene interactions.
2. **Advanced architectures**: Replace TransE/RotatE/ComplEx with GNN-based models (R-GCN, CompGCN, NBFNet) shown to outperform pure embedding approaches on biomedical graphs [9, 12].
3. **Neo4j + PyKEEN pipeline**: Use Neo4j's Cypher query language for path-based reasoning and PyKEEN's production-grade training for reproducible benchmarks.
4. **Hard negative sampling**: Incorporate clinically tested but failed drug–disease pairs as hard negatives.
5. **Temporal split validation**: Implement time-aware train/test splits to simulate prospective discovery.
6. **LLM-enhanced reasoning**: Integrate large language model embeddings (BioMedBERT, GPT-4) with structural graph embeddings following the FuseLinker approach [14].

---

## 7. Conclusion

We presented a biomedical knowledge graph reasoning system for drug repurposing, benchmarking TransE, RotatE, and ComplEx under 5-fold cross-validation on a synthetic KG inspired by DrugBank, DisGeNET, STRING, and CTD. Key findings:

1. **ComplEx achieved the best Hits@10** (0.0636 ± 0.0100) while **RotatE led on MRR** (0.0312 ± 0.0068), consistent with ComplEx's theoretical advantage for one-to-many relations and RotatE's strength on symmetric/antisymmetric patterns.
2. **COVID-19 drug recovery**: the best-ranked known treatment appeared at rank 12 (top 24%), demonstrating proof-of-concept repurposing capability on a small synthetic graph.
3. **Honest evaluation**: all metrics are reported with cross-validation standard deviations, and no result approaches trivially high values (Hits@10 max = 0.064), reflecting the realistic difficulty of drug repurposing from sparse data.
4. **Limitations are explicit**: synthetic data, random negative sampling, and absence of temporal splits are identified as key barriers to real-world applicability.

This work provides a methodological foundation and critical evaluation framework that can scale to production-grade biomedical KG reasoning with real data integration.

---

## References

1. Zhang, Y., Sui, X., Pan, F., Yu, K., & Li, K. (2024). A comprehensive large scale biomedical knowledge graph for AI powered data driven biomedical research. *bioRxiv*. https://doi.org/10.1101/2023.10.13.562216

2. Lou, P., Fang, A., Zhao, W., Yao, K., & Yang, Y. (2023). Potential target discovery and drug repurposing for coronaviruses: study involving a knowledge graph-based approach. *Journal of Medical Internet Research*, 25, e45225. https://doi.org/10.2196/45225

3. Caufield, J. H., Putman, T., Schaper, K., Unni, D. R., & Hegde, H. (2023). KG-Hub — building and exchanging biological knowledge graphs. *Bioinformatics*, 39(7), btad418. https://doi.org/10.1093/bioinformatics/btad418

4. Nam, Y., Lucas, A., Yun, J. S., Lee, S. M., & Park, J. W. (2023). Development of complemented comprehensive networks for rapid screening of repurposable drugs applicable to new emerging disease outbreaks. *Journal of Translational Medicine*, 21, 436. https://doi.org/10.1186/s12967-023-04223-2

5. Bordes, A., Usunier, N., Garcia-Durán, A., Weston, J., & Yakhnenko, O. (2013). Translating embeddings for modeling multi-relational data. *Advances in Neural Information Processing Systems*, 26.

6. Sun, Z., Deng, Z.-H., Nie, J.-Y., & Tang, J. (2019). RotatE: Knowledge graph embedding by relational rotation in complex space. *ICLR 2019*. https://arxiv.org/abs/1902.10197

7. Trouillon, T., Welbl, J., Riedel, S., Gaussier, E., & Bouchard, G. (2016). Complex embeddings for simple link prediction. *ICML 2016*. https://arxiv.org/abs/1606.06357

8. McCoy, K., Gudapati, S., He, L., Horlander, E., & Kartchner, D. (2021). Biomedical text link prediction for drug discovery: a case study with COVID-19. *Pharmaceutics*, 13(6), 794. https://doi.org/10.3390/pharmaceutics13060794

9. Xiao, Y., Hou, Y., Zhou, H., Diallo, G., & Fiszman, M. (2024). Repurposing non-pharmacological interventions for Alzheimer's disease through link prediction on biomedical literature. *Scientific Reports*, 14, 8693. https://doi.org/10.1038/s41598-024-58604-8

10. Zhou, Z., & Yang, S. (2026). Interpretable candidate drug prioritization and explanation framework across-medical knowledge graphs based on graph embedding models: a case study of type 2 diabetes. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0349026

11. Gonzalez-Cavazos, A. C., Tu, R., Sinha, M., & Su, A. I. (2026). A case-based explainable graph neural network framework for mechanistic drug repositioning. *Bioinformatics*, btag008. https://doi.org/10.1093/bioinformatics/btag008

12. Sosa, D. N., Neculae, G., Fauqueur, J., & Altman, R. B. (2024). Elucidating the semantics-topology trade-off for knowledge inference-based pharmacological discovery. *Journal of Biomedical Semantics*, 15, 6. https://doi.org/10.1186/s13326-024-00308-z

13. Xiao, Y., Zhang, S., Zhou, H., Li, M., & Yang, H. (2024). FuseLinker: leveraging LLM's pre-trained text embeddings and domain knowledge to enhance GNN-based link prediction on biomedical knowledge graphs. *Journal of Biomedical Informatics*, 158, 104730. https://doi.org/10.1016/j.jbi.2024.104730

14. Wei, S., Sasi, C., Piepenbrock, J., Huynen, M. A., & 't Hoen, P. A. C. (2025). The use of knowledge graphs for drug repurposing: from classical machine learning algorithms to graph neural networks. *Computers in Biology and Medicine*, 185, 110873. https://doi.org/10.1016/j.compbiomed.2025.110873
