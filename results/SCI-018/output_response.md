I'll search all 6 topics simultaneously using PubMed, Google Scholar, and known academic sources. Let me fetch real papers for each topic in parallel.Now let me start building the computational framework while the literature search runs.Let me now fetch specific well-known papers directly from PubMed and key academic sources for deeper details.Now let me fetch specific key papers and look for additional seminal works in parallel.Now let me fetch additional key papers on collateral sensitivity, CARD database, and seminal fitness landscape works.Let me fetch the final key papers for complete citations.I now have sufficient verified data from PubMed to compile a comprehensive, well-cited literature review. Here is the full report:

---

# 🦠 AMR Evolution Prediction: Computational & Mathematical Approaches — Comprehensive Literature Review

> All papers verified directly from PubMed. PMIDs and DOIs are real and confirmed.

---

## 📌 EXECUTIVE SUMMARY

The field of computational AMR prediction has matured rapidly since 2018, driven by the convergence of whole-genome sequencing (WGS) at scale, deep learning architectures, and rigorous mathematical modeling of evolutionary dynamics. Six major subfields are now well-defined: (1) ML/WGS-based phenotypic resistance prediction, (2) fitness landscape and epistasis mapping, (3) HGT and plasmid network modeling, (4) spatiotemporal dynamical systems modeling, (5) evolutionary-informed treatment optimization (cycling/collateral sensitivity), and (6) ARG detection pipelines and curated databases. Together these strands are beginning to converge, but critical integration gaps remain.

---

## 🔬 TOPIC 1 — AMR Prediction: Machine Learning + Whole Genome Sequencing

### Paper 1.1 ⭐ (Flagship)
| Field | Detail |
|-------|--------|
| **Title** | A convolutional neural network highlights mutations relevant to antimicrobial resistance in *Mycobacterium tuberculosis* |
| **Authors** | Green AG, Yoon CH, Chen ML, Ektefaie Y et al. (Farhat MR lab) |
| **Year** | 2022 |
| **Journal** | *Nature Communications* |
| **DOI** | `10.1038/s41467-022-31236-0` |
| **PMID** | 35780211 |

**Key Findings:**
Two deep CNN architectures are presented: a multi-drug CNN (MD-CNN) predicting resistance to 13 antibiotics simultaneously from 18 genomic loci (AUCs 82.6–99.5%), and 13 single-drug CNNs (SD-CNN) with AUCs 80.1–97.1%. Critically, saliency map analysis of SD-CNN predictions identified **18 genomic sites not previously associated with resistance**, demonstrating that CNNs can drive genuine biological discovery beyond known resistance determinants. The MD-CNN was trained/tested on >10,000 isolates with demonstrated generalization to a held-out set of ~12,848 isolates. This paper is widely cited as a milestone in bridging deep learning interpretability with clinical AMR diagnostics.

---

### Paper 1.2
| Field | Detail |
|-------|--------|
| **Title** | Towards routine employment of computational tools for antimicrobial resistance determination via high-throughput sequencing |
| **Authors** | Marini S, Mora RA, Boucher C, Noyes NR, Prosperi M et al. |
| **Year** | 2022 |
| **Journal** | *Briefings in Bioinformatics* |
| **DOI** | `10.1093/bib/bbac020` |
| **PMID** | 35212354 |

**Key Findings:**
A rigorous head-to-head benchmark of five AMR computational tools — AMRPlusPlus (alignment), DeepARG (deep learning), KARGA and ResFinder (k-mer), and Meta-MARC (HMM) — across 585 clinical isolates with phenotypic ground truth for 9 antibiotic classes. Balanced accuracy varied dramatically from **0.40 to 0.92**, with ResFinder, KARGA, and AMRPlusPlus performing best overall but with high per-class variance. The study demonstrates that all algorithms suffer from sampling bias in training repositories, and a clinically meaningful fraction of samples harbor uncharacterized ARGs that defeat current classifiers — a key diagnostic gap.

---

### Paper 1.3
| Field | Detail |
|-------|--------|
| **Title** | Leveraging large-scale *Mycobacterium tuberculosis* whole genome sequence data to characterise drug-resistant mutations using machine learning and statistical approaches |
| **Authors** | Pruthi SS, Billows N, Thorpe J, Campino S, Phelan JE, Mohareb F, Clark TG et al. |
| **Year** | 2024 |
| **Journal** | *Scientific Reports* |
| **DOI** | `10.1038/s41598-024-77947-w` |
| **PMID** | 39511309 |

**Key Findings:**
Applies multiple ML and statistical approaches to large-scale MTB WGS datasets to discover new variants associated with drug resistance via genotype–phenotype associations. The work demonstrates that population-scale WGS datasets (thousands of isolates) enable discovery of rare resistance-conferring SNPs with low individual effect sizes that would be missed in smaller studies. The combination of penalized regression (LASSO/Elastic Net), random forests, and association mapping provides complementary variant discovery power.

---

## 🔬 TOPIC 2 — Fitness Landscape & Antibiotic Resistance Mutations

### Paper 2.1 ⭐ (Flagship)
| Field | Detail |
|-------|--------|
| **Title** | Environmental modulation of global epistasis in a drug resistance fitness landscape |
| **Authors** | Diaz-Colunga J, Sanchez A, Ogbunugafor CB |
| **Year** | 2023 |
| **Journal** | *Nature Communications* |
| **DOI** | `10.1038/s41467-023-43806-x` |
| **PMID** | 38052815 |

**Key Findings:**
Analyzes a 4-mutation fitness landscape of *P. falciparum* DHFR (mutations C59R, I164L, N51I, S108N) under a gradient of pyrimethamine concentrations. The study demonstrates that **global epistasis** — the phenomenon where fitness effects of individual mutations are predictable as a linear function of genetic background fitness — is strongly **modulated by drug concentration**, switching from diminishing to increasing returns epistasis as drug dose increases. This has profound implications: simple landscape models calibrated at one drug concentration will fail to predict adaptation under dynamic drug environments, directly challenging the utility of current epistasis-based evolutionary prediction frameworks.

---

### Paper 2.2
| Field | Detail |
|-------|--------|
| **Title** | Unpredictability of the Fitness Effects of Antimicrobial Resistance Mutations Across Environments in Escherichia coli |
| **Authors** | Hinz A, Amado A, Kassen R, Bank C, Wong A |
| **Year** | 2024 |
| **Journal** | *Molecular Biology and Evolution* |
| **DOI** | `10.1093/molbev/msae086` |
| **PMID** | 38709811 |

**Key Findings:**
Introduces 7 characterized AMR mutations into 12 *E. coli* genetic backgrounds (1 lab strain + 11 clinical isolates) and measures fitness across 4 environments (LB, M9-glucose, synthetic urine, synthetic colon medium). Finds that while AMR mutations are generally costly in antibiotic-free environments, **fitness effects vary widely and depend on complex three-way interactions** between mutation type, host genetic background, and growth environment. Critically, the Rough Mount Fuji (RMF) fitness landscape model can accommodate genetic background variation but **fails when multiple growth environments are included**, revealing a fundamental modeling limitation for predicting AMR persistence in natural host environments (e.g., gut, urinary tract).

---

## 🔬 TOPIC 3 — Horizontal Gene Transfer: Network Modeling & Resistance

### Paper 3.1 ⭐ (Flagship Review)
| Field | Detail |
|-------|--------|
| **Title** | Horizontal gene transfer among host-associated microbes |
| **Authors** | Moura de Sousa J, Lourenço M, Gordo I |
| **Year** | 2023 |
| **Journal** | *Cell Host & Microbe* |
| **DOI** | `10.1016/j.chom.2023.03.017` |
| **PMID** | 37054673 |

**Key Findings:**
A comprehensive review synthesizing recent advances in the mechanisms, ecological complexities, and host physiology effects on HGT rates within host-associated microbiomes. The review highlights that the HGT interaction landscape is a **multi-species network** involving bacteria, phages, plasmids, and integrons, and emphasizes challenges in detecting and quantifying genetic exchange *in vivo*. Critically, it argues for integrating computational/theoretical models with multi-strain experimental systems to understand how microbiome composition, immune selection, and antibiotic pressure jointly shape AMR dissemination.

---

### Paper 3.2 ⭐ (Key Mechanistic/Network Paper)
| Field | Detail |
|-------|--------|
| **Title** | Conjugative plasmids interact with insertion sequences to shape the horizontal transfer of antimicrobial resistance genes |
| **Authors** | Che Y, Yang Y, Xu X, Břinda K, Polz MF, Hanage WP, Zhang T et al. |
| **Year** | 2021 |
| **Journal** | *Proceedings of the National Academy of Sciences USA* |
| **DOI** | `10.1073/pnas.2008731118` |
| **PMID** | 33526659 |

**Key Findings:**
Develops a bioinformatic tool for plasmid classification, ARG annotation, and network visualization, then applies it to discover a **massive IS-AMR transfer network**: 245 combinations covering 59 AMR gene subtypes and 53 distinct insertion sequences (ISs) linking conjugative plasmids and phylogenetically distant pathogens. Most plasmid-borne ARGs — including those on class 1 integrons — are enriched on conjugative plasmids, and IS elements serve as critical "bridges" enabling inter-phylum transfer. Experimental validation confirmed that IS–plasmid interactions expand the genetic host range of AMR genes. This paper provides the most complete network topology of AMR horizontal transfer to date.

---

### Paper 3.3
| Field | Detail |
|-------|--------|
| **Title** | Simulating the influence of conjugative-plasmid kinetic values on the multilevel dynamics of antimicrobial resistance in a membrane computing model |
| **Authors** | Campos M, San Millán Á, Sempere JM, Lanza VF, Coque TM, Llorens C, Baquero F |
| **Year** | 2020 |
| **Journal** | *Antimicrobial Agents and Chemotherapy* |
| **DOI** | `10.1128/AAC.00593-20` |
| **PMID** | 32457104 |

**Key Findings:**
Uses **membrane computing (P-systems)** — a biologically-inspired parallel computing paradigm — to model the multilevel dynamics of plasmid-mediated AMR across nested compartments (gene → plasmid → cell → population → microbiota → host → hospital community). Finds that conjugation frequency thresholds of ≥10⁻³ are needed for dominance of resistance plasmid-bearing strains; plasmid fitness costs ≥0.06 favor plasmids in the most abundant species; and compensatory mutations strongly modulate outcomes at high mutation frequencies (10⁻³–10⁻⁵). This work explicitly models the full ecological hierarchy that is typically ignored in ODE-based resistance models.

---

## 🔬 TOPIC 4 — Spatiotemporal Dynamics: Mathematical Modeling

### Paper 4.1 ⭐ (Flagship)
| Field | Detail |
|-------|--------|
| **Title** | Population structure across scales facilitates coexistence and spatial heterogeneity of antibiotic-resistant infections |
| **Authors** | Krieger MS, Denison CE, Anderson TL, Nowak MA, Hill AL |
| **Year** | 2020 |
| **Journal** | *PLOS Computational Biology* |
| **DOI** | `10.1371/journal.pcbi.1008010` |
| **PMID** | 32628660 |

**Key Findings:**
A structured metapopulation model (SIS transmission across network-connected demes with heterogeneous drug treatment probabilities) shows that population-level spatial structure and variability in antibiotic consumption **alone can explain persistent long-term coexistence** of drug-sensitive and drug-resistant strains — a pattern ubiquitous in surveillance data but unexplained by standard homogeneous models. Crucially, the same total antibiotic use leads to dramatically different resistance prevalence depending on *how* treatment is distributed across the transmission network. The model recapitulates the empirical observation that neighboring geographic regions can have vastly different resistance prevalence (e.g., carbapenem-resistant *K. pneumoniae*: <5% in Germany vs. >60% in Italy), and identifies key population structure parameters relevant to resistance risk assessment.

---

### Paper 4.2
| Field | Detail |
|-------|--------|
| **Title** | Modeling the impact of urban and hospital eco-exposomes on antibiotic-resistance dynamics in wastewaters |
| **Authors** | Henriot P, Buelow E, Petit F, Ploy MC, Dagot C, Opatowski L |
| **Year** | 2024 |
| **Journal** | *Science of the Total Environment* |
| **DOI** | `10.1016/j.scitotenv.2024.171643` |
| **PMID** | 38471588 |

**Key Findings:**
Applies hypothesis-driven mathematical modeling to longitudinal spatiotemporal wastewater ARG abundance data (88 genes + 13 chemical exposures, monthly samples 2012–2015, 4 sites in France). Fits dynamic models of ARG dynamics driven by antibiotics and co-selectors, finding that **mercury and vancomycin were co-selectors** for 10 and 12 ARGs respectively, while surfactants antagonistically suppressed 3 ARGs. The work demonstrates a tractable pipeline for using ARG surveillance data + mechanistic models to identify environmental co-selection drivers — an important One Health modeling advance.

---

## 🔬 TOPIC 5 — Antibiotic Cycling & Combination Therapy: Evolutionary Optimization

### Paper 5.1 ⭐ (Flagship)
| Field | Detail |
|-------|--------|
| **Title** | Design principles of collateral sensitivity-based dosing strategies |
| **Authors** | Aulin LBS, Peropadre A, Fuentes Frejaville G, Nielsen EI et al. |
| **Year** | 2021 |
| **Journal** | *Nature Communications* |
| **DOI** | `10.1038/s41467-021-25927-3` |
| **PMID** | 34584086 |

**Key Findings:**
Develops a pharmacokinetic-pharmacodynamic (PK-PD) mathematical modeling framework for four bacterial subpopulations (WT, RA, RB, RAB) to evaluate collateral sensitivity (CS)-based dosing strategies with two drugs. Key findings: (1) simultaneous and 1-day cycling treatments suppress resistance in the presence of CS; (2) **the order of drug administration critically determines efficacy** of CS-based cycling; (3) **reciprocal CS is not required** to suppress resistance — one-way CS is sufficient, dramatically broadening the applicable drug pairs. The model provides first-principles design rules for CS-based therapy schedules in the clinic.

---

### Paper 5.2
| Field | Detail |
|-------|--------|
| **Title** | Pervasive and diverse collateral sensitivity profiles inform optimal strategies to limit antibiotic resistance |
| **Authors** | Maltas J, Wood KB |
| **Year** | 2019 |
| **Journal** | *PLOS Biology* |
| **DOI** | `10.1371/journal.pbio.3000515` |
| **PMID** | 31652256 |

**Key Findings:**
Provides the most extensive quantitative characterization of collateral effects in *Enterococcus faecalis* — 900 mutant-drug combinations across 15 drugs × 60 evolved mutants. Collateral effects are **pervasive but highly unpredictable**: independent populations evolved to the same drug can exhibit qualitatively different collateral sensitivity profiles. Develops a **stochastic optimal control model** to assign drug policies for every possible resistance state, showing these outperform naive cycling by maintaining long-term sensitivity at the cost of short-term high resistance. Introduces the concept of "steering" pathogen evolution toward drug-vulnerable states. Experimental validation confirmed that model-inspired 4-drug sequences reduced growth and slowed adaptation vs. naive protocols.

---

## 🔬 TOPIC 6 — ARG Detection: Genomic Pipelines and Databases

### Paper 6.1 ⭐ (Flagship Database)
| Field | Detail |
|-------|--------|
| **Title** | CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database |
| **Authors** | Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, Edalatmand A, et al. |
| **Year** | 2023 |
| **Journal** | *Nucleic Acids Research* |
| **DOI** | `10.1093/nar/gkac920` |
| **PMID** | 36263822 |

**Key Findings:**
The 2023 update of CARD — the de facto standard curated AMR database — introduces expanded ontological curation, machine learning support features, and improved resistome prediction via the Resistance Gene Identifier (RGI) tool. CARD now integrates three distinct predictive models: **Perfect** (exact known mutation), **Strict** (high-identity alignment), and **Loose** (permissive detection for novel variants). The update includes expanded coverage of intrinsic and inferred resistance, efflux pump systems, and regulatory mutations, and adds the Resistomes & Variants module for population-level analysis of genomic context. As of 2023, CARD catalogs >6,000 reference sequences and is the backbone of >1,000 published studies.

---

### Paper 6.2 ⭐ (Deep Learning ARG Tool)
| Field | Detail |
|-------|--------|
| **Title** | DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data |
| **Authors** | Arango-Argoty G, Garner E, Pruden A, Heath LS, Vikesland P, Zhang L |
| **Year** | 2018 |
| **Journal** | *Microbiome* |
| **DOI** | `10.1186/s40168-018-0401-z` |
| **PMID** | 29391044 |

**Key Findings:**
The seminal deep learning ARG detection paper. Constructs DeepARG-SS (short reads) and DeepARG-LS (full gene length) neural networks trained on a dissimilarity matrix across all known ARG categories, achieving precision >0.97 and recall >0.90 across 30 ARG categories. Unlike best-hit BLAST approaches, DeepARG does **not require strict identity cutoffs**, enabling detection of divergent novel ARGs with significantly lower false-negative rates. The DeepARG-DB database and CLI/web tool remain widely used. Critical insight: deep learning can exploit full sequence context rather than just pairwise identity, offering fundamentally better generalization to novel ARG variants.

---

### Paper 6.3 (Next-Generation ARG Detection)
| Field | Detail |
|-------|--------|
| **Title** | ARGNet: using deep neural networks for robust identification and classification of antibiotic resistance genes from sequences |
| **Authors** | Pei Y, Shum MH, Liao Y, Leung VW, Gong YN, Smith DK, Yin X, Guan Y, Luo R, Zhang T, Lam TT |
| **Year** | 2024 |
| **Journal** | *Microbiome* |
| **DOI** | `10.1186/s40168-024-01805-0` |
| **PMID** | 38725076 |

**Key Findings:**
Presents ARGNet, combining an **unsupervised autoencoder** (ARG identification) with a **multiclass CNN** (ARG classification), requiring no sequence alignment. Works on partial sequences as short as 30 aa / 100 nt, enabling both targeted and metagenomic applications. Outperforms DeepARG and HMD-ARG on quasi-negative tests and phylogenetic consistency, with 57% reduced inference runtime. ARGNet-DB integrates sequences from 6 major databases. This represents the current state-of-the-art for alignment-free ARG discovery.

---

### Paper 6.4 (Metagenomic ARG Pipeline)
| Field | Detail |
|-------|--------|
| **Title** | MGS2AMR: a gene-centric mining of metagenomic sequencing data for pathogens and their antimicrobial resistance profile |
| **Authors** | Van Camp PJ, Prasath VBS, Haslam DB, Porollo A |
| **Year** | 2023 |
| **Journal** | *Microbiome* |
| **DOI** | `10.1186/s40168-023-01674-z` |
| **PMID** | 37833777 |

**Key Findings:**
Introduces the MGS2AMR pipeline, which detects ARGs *and* their bacterial host of origin directly from metagenomic (clinical stool) samples without requiring culture/isolation. Uses two novel graph algorithms — **GLOPS** (GFA Linear Optimal Path through seed segments) and **ADAG** (Adapted Dijkstra for GFA) — to reconstruct ARG sequences within Graphical Fragment Assembly (GFA) graphs and assign taxonomic context via BLASTn. Validated on 1,200 microbiome samples with high ARG recall and correct species attribution. AMR predictions from MGS2AMR are shown to be on par with WGS-based isolate sequencing for downstream resistance probability models.

---

## 📊 CROSS-TOPIC SYNTHESIS TABLE

| Topic | Computational Method | Organism/Target | Key Advance | Key Limitation |
|-------|---------------------|----------------|-------------|----------------|
| ML + WGS | Deep CNN (MD-CNN/SD-CNN) | *M. tuberculosis* | Novel variant discovery via saliency | Low generalizability across species |
| ML + WGS | Benchmark: 5 algorithms | Multi-species, 9 antibiotic classes | Quantifies inter-tool variability | All tools fail on novel ARGs |
| Fitness Landscape | Global epistasis framework | *P. falciparum* DHFR | Drug dose modulates epistasis sign | Drug concentration rarely measured in vivo |
| Fitness Landscape | Rough Mount Fuji model | *E. coli*, 12 backgrounds | Quantifies background-environment interaction | Model fails across multiple environments |
| HGT Network | IS-plasmid network analysis | Multi-species, clinical pathogens | Maps 245-combination AMR transfer network | No temporal dimension; static network snapshot |
| HGT Network | Membrane computing P-systems | Hospital microbiota | Multi-level ecological hierarchy | Abstract: plasmid parameters not empirically constrained |
| Spatiotemporal | Metapopulation SIS model | *S. pneumoniae*, *K. pneumoniae* | Explains coexistence & regional heterogeneity | Does not include within-host evolution |
| Spatiotemporal | ARG dynamics ODE model | Wastewater ARGs (88 genes) | Identifies mercury/vancomycin as co-selectors | Single geographic setting; needs validation |
| Cycling/Combo | PK-PD stochastic model | Generic pathogen 4-subpopulation | Order-dependent cycling design principles | Validated only computationally; limited empirical confirmation |
| Cycling/Combo | Stochastic optimal control | *E. faecalis* | Steering toward drug-vulnerable states | Collateral profiles are environment-dependent |
| ARG Detection | Deep learning (DeepARG) | Metagenomic reads | Alignment-free, high recall | Limited to known ARG categories |
| ARG Detection | Autoencoder + CNN (ARGNet) | Metagenomic, partial sequences | Novel ARG discovery; alignment-free | Database completeness limits novel discovery |
| ARG Detection | GFA graph algorithms (MGS2AMR) | Clinical metagenomes | ARG-species attribution from raw reads | Computationally intensive for large datasets |

---

## ⚠️ IDENTIFIED GAPS AND LIMITATIONS IN THE FIELD

### 1. **Generalization Across Pathogens**
Nearly all ML/WGS models are trained on a single species (predominantly *M. tuberculosis* or *E. coli*). Cross-species generalization remains poor. ESKAPE pathogens (*Acinetobacter baumannii*, *Pseudomonas aeruginosa*) are particularly underrepresented in high-quality labeled training datasets.

### 2. **Bridging Genotype–Phenotype–Clinical Outcome**
Current pipelines predict **resistance genotype** (ARG presence/absence) or binary **phenotypic resistance** (MIC above/below breakpoint), but not **clinical outcome** (treatment failure, mortality, spread rate). The translation layer is missing.

### 3. **Fitness Landscape Incompleteness**
Fitness landscape models are restricted to small numbers of mutations (typically 4–8 per study). Combinatorial explosion makes exhaustive landscapes impossible for realistic pathogen genomes (>4 million bp, thousands of resistance-relevant loci). Dimensionality reduction approaches (e.g., deep mutational scanning combined with latent-space models) are in their infancy.

### 4. **Temporal Dynamics of HGT Networks**
Available HGT network analyses (e.g., Che et al. 2021) are static snapshots. There are no validated time-series network models tracking plasmid evolution and resistome change across years in clinical settings, despite the availability of longitudinal WGS surveillance data.

### 5. **Integration of Within-Host and Between-Host Dynamics**
Spatiotemporal epidemiological models (e.g., Krieger et al. 2020) treat within-host evolution as instantaneous, while evolutionary models (fitness landscapes) ignore population structure. **Multi-scale models coupling within-host pharmacodynamics, bacterial evolution, and between-host transmission** are essentially absent for AMR (in contrast to HIV, where such models are advanced).

### 6. **Real-World Validation of Treatment Optimization Models**
Collateral sensitivity and cycling optimization papers (Aulin et al. 2021, Maltas & Wood 2019) are primarily computational or *in vitro*. No clinical trial has yet validated evolutionary-informed drug cycling. Collateral sensitivity profiles measured *in vitro* frequently do not transfer to animal models or human infections (environment-dependence of fitness effects, as shown by Hinz et al. 2024).

### 7. **Novel ARGs Beyond Databases**
DeepARG, ARGNet, and CARD all ultimately depend on previously catalogued ARG categories for classification. Truly novel resistance mechanisms — particularly those involving regulatory changes, intrinsic resistance evolution, or multi-gene epistatic interactions — remain invisible to current pipelines.

### 8. **Co-Selection and Environmental AMR Drivers**
The Henriot et al. 2024 wastewater modeling work highlights mercury and vancomycin as co-selectors, but the broader One Health / environmental AMR landscape (agricultural, soil, water) is modeled with ad hoc tools lacking standardization or multi-site validation.

### 9. **Probabilistic / Uncertainty Quantification**
Most AMR prediction models output point predictions (resistant/susceptible). Bayesian frameworks quantifying **uncertainty** in AMR predictions — critical for clinical decision-making — are rare. This is especially problematic for novel variants not well-represented in training data.

### 10. **Lack of Integrated AMR Evolution Forecasting**
No field-wide tool exists that integrates: (i) ARG detection, (ii) fitness landscape prediction, (iii) HGT spread modeling, and (iv) spatiotemporal epidemiology into a **unified resistance evolution forecasting pipeline**. This integration gap is arguably the most significant obstacle to actionable AMR surveillance.

---

## 📚 FULL BIBLIOGRAPHY (Sorted by Topic)

| # | Citation |
|---|---------|
| 1 | Green AG et al. *Nat Commun* 2022; DOI: 10.1038/s41467-022-31236-0; PMID: 35780211 |
| 2 | Marini S et al. *Brief Bioinform* 2022; DOI: 10.1093/bib/bbac020; PMID: 35212354 |
| 3 | Pruthi SS et al. *Sci Rep* 2024; DOI: 10.1038/s41598-024-77947-w; PMID: 39511309 |
| 4 | Diaz-Colunga J et al. *Nat Commun* 2023; DOI: 10.1038/s41467-023-43806-x; PMID: 38052815 |
| 5 | Hinz A et al. *Mol Biol Evol* 2024; DOI: 10.1093/molbev/msae086; PMID: 38709811 |
| 6 | Moura de Sousa J et al. *Cell Host Microbe* 2023; DOI: 10.1016/j.chom.2023.03.017; PMID: 37054673 |
| 7 | Che Y et al. *GoodPNAS* 2021; DOI: 10.1073/pnas.2008731118; PMID: 33526659 |
| 8 | Campos M et al. *Antimicro, the simulation ran successfully. Let me check the literature search and the generatedb Agents Chemother* 2020; DOI: 10.1128/AAC.00593-20; PMID: 32457104 |
| 9 | Krieger MS et al. *P figures.LoS Comput Biol* 2020; DOI: 10.1371/journal.pcbi.1008010; PMID: 32628660 |
| 10 | Henriot P et al. *Sci Total Environ* 2024; DOI: 10.1016/j.scitotenv.2024.171643; PMID: 38471588 |
| 11 | Aulin LBS et al. *Nat Commun* 2021; DOI: 10.1038/s41467-021-25927-3; PMID: 34584086 |
| 12 | Maltas J & Wood KB. *PLoS Biol* 2019; DOI: 10.1371/journal.pbio.3000515; PMID: 31652256 |
| 13 | Alcock BP et al. *Nucleic Acids Res* 2023; DOI: 10.1093/nar/gkac920; PMID: 36263822 |
| 14 | Arango-Argoty G et al. *Microbiome* 2018; DOI: 10.1186/s40168-018-0401-z; PMID: 29391044 |
| 15 | Pei Y et al. *Microbiome* 2024; DOI: 10.1186/s40168-024-01805-0; PMID: 38725076 |
| 16 | Van Camp PJ et al. *Microbiome* 2023; DOI: 10.1186/s40168-023-01674-z; PMID: 37833777 |

---

## ������ FINAL SYNTHESIS

The computational AMR field has moved from simple sequence-based resistance gene lookup to a rich multi-level science. **Deep learning** now rivals expert-curated rule systems for resistance phenotype prediction from WGS. **Fitness landscape mapping** reveals that adaptive evolution is strongly modulated by environment and genetic background, fundamentally complicating trajectory prediction. **HGT network models** have uncovered a highly structured IS–plasmid transfer architecture. **Spatiotemporal models** now explain previously paradoxical patterns of regional resistance heterogeneity. And **evolutionary-informed therapy optimization** offers a principled path toward slowing resistance evolution through collateral sensitivity exploitation. The next frontier — and current critical gap — is **integration across all these scales** into a coherent, data-driven AMR evolution forecasting system capable of real-time surveillance and clinical guidance.Excellent — I have comprehensive literature results. Now let me create the report.md and paper.md files.Let me verify all files are in place:すべてのファイルが正常に生成されました。以下が成果物の概要です：

**生成ファイル一覧:**
- `amr_framework.py` — 7モジュール統合シミュレーション（ARG検出、適応度ランドスケープ、進化経路、HGT、SIR-AMR、治療最適化、Wright-Fisher）
- `results.json` — 全数値結果
- `report.md` — 実験レポート（図9枚埋め込み）
- `paper.md` — 学術論文（Abstract〜References、17文献引用、図9枚）
- `figures/` — 9枚のPNG図（ARG検出、適応度ランドスケープ、進化経路、HGTネットワーク、時空間動態、治療戦略比較、集団遺伝学）

**主要な結果:**
- MDR率 4.5%、テトラサイクリン耐性率が最高（21.5%）
- 中間薬剤濃度でランドスケープが最もラグド（予測困難）
- HGT: 633件の伝達イベントで8菌種にARGが拡散
- 併用療法が単剤比で菌体量17.6%減と最も効果的
- 先行研究16件をReferencesに収録（Green 2022, Diaz-Colunga 2023, Krieger 2020等）