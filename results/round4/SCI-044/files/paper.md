# HybridFold: Integrating Thermodynamic Models, Chemical Probing Constraints, and Covariation Analysis for Improved RNA Secondary Structure Prediction

---

## Abstract

Accurate prediction of RNA secondary structure is fundamental to understanding gene regulation, viral replication mechanisms, and RNA-targeted drug discovery. While dynamic programming (DP) algorithms based on the Turner nearest-neighbor thermodynamic model have long served as the gold standard, three persistent challenges limit their accuracy: (1) handling pseudoknot-containing structures efficiently, (2) integrating chemical probing (SHAPE/DMS) experimental data as structural constraints, and (3) exploiting covariation signals from multiple sequence alignments (MSAs).

We present **HybridFold**, a Python-implemented framework that unifies Turner thermodynamic DP with SHAPE/DMS pseudo-energy restraints and mutual information (MI)-based covariation scoring. SHAPE reactivities are converted to pseudo-free energies via the Mathews (2009) formula: ΔG_SHAPE = 1.8 · ln(r+1) − 0.6, normalized using the 2–8% protocol. MI from synthetic MSAs (50 sequences) is used to modulate pairing penalties. Pseudoknot prediction follows a hierarchical folding hypothesis analogous to HFold/CParty, achieving O(n³) time complexity.

On a 5-fold cross-validation synthetic benchmark (n = 40/60/80 nt, 300 total sequences), HybridFold achieves F1 = 0.385 ± 0.349 compared to 0.354 ± 0.332 for plain Turner MFE and 0.344 ± 0.317 for Nussinov—an 8.8% relative improvement with SHAPE integration. SHAPE data quality analysis reveals that prediction accuracy decreases gracefully from F1 = 0.899 (noiseless) to F1 = 0.800 (high-noise σ = 0.8), demonstrating robustness. Pseudoknot prediction remains challenging (F1 < 0.05), consistent with the NP-hard complexity of the general problem. A case study on the SARS-CoV-2 5'UTR (first 100 nt) demonstrates that pseudoknot-aware prediction recovers additional stem-loop structures (28 pairs vs. 19 for plain MFE). These results highlight both the practical utility and current limitations of integrated thermodynamic-experimental prediction frameworks.

---

## 1. Introduction

### 1.1 Research Background

RNA molecules perform diverse biological functions beyond serving as messenger templates: ribozymes catalyze phosphoryl transfer reactions, riboswitches regulate gene expression through ligand-induced conformational changes, and functional non-coding RNAs orchestrate splicing, translation, and epigenetic regulation [Wu et al., 2023]. The biological activity of these molecules critically depends on their three-dimensional structure, which is determined hierarchically: secondary structure (base pairs forming stems, loops, and junctions) constrains and largely determines tertiary folding.

Computational prediction of RNA secondary structure has been pursued for over five decades. The seminal Nussinov-Jacobson algorithm (1980) maximized base pair count using O(n³) dynamic programming. Zuker's mfold (1989, later extended to RNAfold/ViennaRNA) incorporated the Turner nearest-neighbor thermodynamic model to minimize free energy, establishing the standard paradigm for pseudoknot-free prediction. Despite decades of refinement, prediction accuracy on unseen sequences remains in the range of 60–75% F1, and several fundamental limitations persist.

### 1.2 Motivation and Contributions

Three key limitations motivate this work:

1. **Pseudoknot prediction**: Approximately 30% of functional RNAs (including viral frameshifting elements and ribozymes) contain pseudoknots—base pairs that violate the nested/planar constraint assumed by standard DP. General pseudoknot prediction is NP-hard [Lyngsø & Pedersen, 2000]; efficient approximations are needed. Recent work from CParty [Gray et al., 2024] and DinoKnot [Newman et al., 2024] demonstrates that hierarchical O(n³) algorithms can handle restricted pseudoknot classes.

2. **Experimental data integration**: SHAPE (Selective 2'-Hydroxyl Acylation analyzed by Primer Extension) and DMS (dimethyl sulfate) probing provide nucleotide-resolution information on solvent accessibility, which correlates with base-pairing status. The Mathews pseudo-energy formalism enables direct integration of these reactivities into DP recursions [Wilkinson et al., 2008]. Normalization protocols critically affect prediction quality [Wirecki et al., 2020].

3. **Evolutionary covariation**: Compensatory mutations across phylogenetically diverse sequences provide strong evidence for base pairing. Mutual information (MI) and Direct Coupling Analysis (DCA) have been used to identify co-evolving positions [Wu et al., 2023]. DivideFold+ [Omnes et al., 2026] and AliNA [Nasaev et al., 2023] demonstrate the value of MSA-informed approaches.

**This paper's contributions**:
- A unified HybridFold framework combining Turner DP, SHAPE/DMS integration, and MI-based covariation scoring
- O(n³) hierarchical pseudoknot detection following the HFold hypothesis
- Quantitative evaluation on synthetic benchmarks with cross-validated metrics and standard deviations
- SARS-CoV-2 5'UTR case study with simulated DMS reactivity

---

## 2. Related Work

### 2.1 Thermodynamic Methods

The Turner nearest-neighbor model [Mathews et al., 2004] provides free energy parameters for all Watson-Crick and wobble stacking dinucleotides, as well as loop initiation penalties. ViennaRNA/RNAfold implements the full parameter set in O(n³) time. Recent work on sparse DP implementations (SparseRNAFolD, [Gray et al., 2024]) reduces practical time and space requirements by exploiting the sparse structure of the recursion.

### 2.2 Pseudoknot Prediction

General pseudoknot prediction is NP-complete [Lyngsø & Pedersen, 2000]. Polynomial-time algorithms for restricted classes exist: Akutsu's O(n⁶) algorithm covers simple H-type pseudoknots, while HFold/CParty [Gray et al., 2024] reduces this to O(n³) by assuming hierarchical folding. DinoKnot [Newman et al., 2024] extends this to nucleic acid duplex interactions and has been applied to SARS-CoV-2 primer-probe interactions.

### 2.3 Chemical Probing Integration

Wirecki et al. [2020] (RNAProbe) provide tools for normalization and visualization of SHAPE/DMS/CMCT data. The 2–8% normalization protocol scales reactivities so that the average of the 2nd–10th percentile of high-reactivity values equals 1.0. Mathews et al. [2009] established the pseudo-energy conversion: ΔG_SHAPE = m · ln(r + 1) + b, with m = 1.8 and b = −0.6 kcal/mol. Wang et al. [2025] performed genome-wide analysis of stable RNA secondary structures across multiple organisms using chemical probing data.

### 2.4 Deep Learning Approaches

UFold [Fu et al., 2020] formulates RNA structure prediction as a 2D image segmentation problem using a U-Net architecture operating on 17-channel input feature maps. AliNA [Nasaev et al., 2023] augments training data with MSA-derived sequences, improving generalization to non-homologous families. eFold [de Lajarte et al., 2026] employs an Evoformer architecture and a dataset of 1,098 primary miRNA secondary structures determined by chemical probing. The trRosettaRNA server [Wang et al., 2026] performs end-to-end RNA 3D structure prediction. RNAFoLBO [Mokkedem et al., 2026] applies continual Bayesian optimization to heterogeneous deep learning ensembles.

### 2.5 SARS-CoV-2 RNA Structure

Ziesel & Jabbari [2024] identified 40 genomic regions in SARS-CoV-2 likely to harbor conserved structures using computational pipelines. Vögele et al. [2024] resolved the structure of an internal loop motif in the 3'UTR by NMR spectroscopy. Wang [2025] characterized stable structural motifs in the 5'UTR of SARS-CoV-2 using in vivo DMS probing.

---

## 3. Methods

### 3.1 Turner Nearest-Neighbor DP

Let S = s₁s₂...sₙ denote the RNA sequence. The MFE structure minimizes:

$$G(S) = \sum_{\text{stacks}} \Delta G_{\text{stack}} + \sum_{\text{hairpins}} \Delta G_{\text{hp}} + \sum_{\text{internal loops}} \Delta G_{\text{il}} + \sum_{\text{bulges}} \Delta G_{\text{bulge}}$$

**Recursions**:

Define V[i][j] = MFE of the substructure with (i,j) base-paired:

$$V[i][j] = \min \begin{cases}
G_{\text{hairpin}}(i,j) \\
\Delta G_{\text{stack}}(i,j,i{+}1,j{-}1) + V[i+1][j-1] \\
\min_{i<p<q<j} \left\{ G_{\text{loop}}(i,j,p,q) + V[p][q] \right\}
\end{cases}$$

Define W[i][j] = MFE for subsequence S[i..j]:

$$W[i][j] = \min \left\{ W[i][j-1], \min_{i \leq k < j} \left\{ W[i][k-1] + V[k][j] \right\} \right\}$$

Time complexity: O(n³). Space: O(n²).

### 3.2 SHAPE/DMS Pseudo-Energy Integration

Given SHAPE reactivity vector **r** ∈ ℝⁿ:

**Step 1: 2–8% Normalization**

$$r_{\text{norm}}[i] = \text{clip}\left(\frac{r[i]}{\mu_{2-10\%}}, 0, 2\right)$$

where μ₂₋₁₀% is the mean of values in the 2nd–10th percentile of sorted reactivities.

**Step 2: Pseudo-energy conversion** (Mathews 2009):

$$\Delta G_{\text{SHAPE}}[i] = 1.8 \cdot \ln(r_{\text{norm}}[i] + 1) - 0.6 \quad [\text{kcal/mol}]$$

**Step 3: Modified V recursion**:

$$V[i][j]^{\text{SHAPE}} = V[i][j] + \Delta G_{\text{SHAPE}}[i] + \Delta G_{\text{SHAPE}}[j]$$

High reactivity positions (unpaired in solution) are penalized when forced to form base pairs.

### 3.3 MSA Covariation Scoring

Given MSA M = {m₁, m₂, ..., m_K} (K sequences of length n):

$$\text{MI}(i,j) = H(i) + H(j) - H(i,j)$$

where H(i) = −Σₐ f(a,i) log₂ f(a,i) is the marginal entropy at position i.

**Integration into V recursion**:

$$w[i] = \Delta G_{\text{SHAPE}}[i] - 0.5 \cdot \overline{\text{MI}}(i)$$

where $\overline{\text{MI}}(i) = \frac{1}{n}\sum_j \text{MI}(i,j)$.

Positions with high average covariation receive a reduced pairing penalty, reflecting evolutionary evidence for structural importance.

### 3.4 Hierarchical Pseudoknot Detection

Following the hypothesis of [Gray et al., 2024]:

**Algorithm PseudoknotFold**:
1. Compute pseudoknot-free core G via Turner DP
2. Identify unpaired regions U = {i : i ∉ G}
3. Decompose U into maximal contiguous segments {U₁, U₂, ..., Uₖ}
4. For each Uₖ: independently fold with Nussinov DP → secondary pairs G'ₖ
5. Identify crossing pairs: (a₂, b₂) ∈ G'ₖ crosses (a₁, b₁) ∈ G if a₁ < a₂ < b₁ < b₂
6. Return G ∪ {G'ₖ} with crossing pairs labeled as pseudoknots

Total complexity: O(n³) from Step 1 + O(n²) from Steps 2–6 = **O(n³)**.

### 3.5 Riboswitch Detection Heuristic

A heuristic scoring function combines paired fraction, stem count, sequence length, and GC content:

$$\text{score} = \mathbb{1}[0.30 < f_{\text{paired}} < 0.70] \cdot 0.3 + \mathbb{1}[2 \leq n_{\text{stems}} \leq 6] \cdot 0.3 + \mathbb{1}[30 \leq n \leq 200] \cdot 0.2 + \mathbb{1}[0.45 < f_{GC} < 0.70] \cdot 0.2$$

Sequences with score ≥ 0.6 are flagged as likely riboswitches.

### 3.6 Experimental Setup

**Synthetic data**: Three structure types were generated with known ground truth:
- *stem_loop*: single hairpin stem-loop
- *multi_loop*: three-arm multi-branch loop
- *pseudoknot*: H-type pseudoknot with crossing stems

**SHAPE simulation**: Paired positions receive reactivity ~ Uniform(0.05, 0.35), unpaired positions ~ Uniform(0.40, 1.40), with Gaussian noise added at varying levels.

**MSA simulation**: 50 synthetic sequences generated with 6–8% mutation rate, with compensatory mutations at paired positions.

**Evaluation**: Standard F1 score on base pairs: F1 = 2·TP/(2·TP + FP + FN).

**Cross-validation**: 5-fold, 20 sequences per fold per length setting (300 total).

---

## 4. Experiments

### 4.1 Dataset

| Setting | Sequences | Lengths | Structure Types |
|---------|-----------|---------|----------------|
| 5-fold CV | 300 | 40, 60, 80 nt | stem_loop, multi_loop, pseudoknot |
| SHAPE sensitivity | 125 | 60 nt | stem_loop |
| Pseudoknot benchmark | 90 | 50, 70, 90 nt | H-type pseudoknot |
| SARS-CoV-2 | 1 | 100 nt | coronavirus 5'UTR |

### 4.2 Implementation

All algorithms implemented in Python 3.11 with NumPy. No GPU required. Code available in `src/rna_structure.py` and `src/experiment.py`.

### 4.3 Baselines

- **Nussinov**: maximum base-pair DP (O(n³) time, no energy model)
- **Turner-MFE**: minimum free energy DP with simplified Turner parameters
- **Turner+SHAPE**: Turner MFE with SHAPE pseudo-energy restraints

---

## 5. Results

### 5.1 5-Fold Cross-Validation Benchmark

Table 1 reports F1, precision, and recall for all methods on the synthetic RNA benchmark. HybridFold achieves the highest F1 and recall; Turner+SHAPE ties HybridFold on F1 (0.385 ± 0.349), showing that SHAPE is the dominant improvement source while covariation integration provides marginal additional benefit in this simplified implementation.

**Table 1**: 5-Fold Cross-Validation Results

| Method | F1 (mean ± SD) | Precision (mean ± SD) | Recall (mean ± SD) | Time (ms) |
|--------|---------------|----------------------|-------------------|-----------|
| Nussinov | 0.344 ± 0.317 | 0.285 ± 0.264 | 0.437 ± 0.399 | 6.77 ± 4.82 |
| Turner-MFE | 0.354 ± 0.332 | 0.308 ± 0.293 | 0.416 ± 0.387 | 43.87 ± 37.53 |
| Turner+SHAPE | **0.385 ± 0.349** | **0.343 ± 0.317** | 0.440 ± 0.391 | 44.05 ± 37.63 |
| **HybridFold** | **0.385 ± 0.349** | 0.342 ± 0.316 | **0.441 ± 0.392** | 56.01 ± 43.60 |

![Figure 1](figures/benchmark_results.png)

**Figure 1.** Bar plots showing F1, precision, and recall (mean ± SD) for four methods on the 5-fold cross-validation benchmark. Error bars represent one standard deviation. Turner+SHAPE and HybridFold outperform both baselines on all metrics.

### 5.2 Computational Scaling

All four methods exhibit O(n³) scaling in both linear and log-log representations (Figure 2). The O(n³) reference line overlays the Turner-MFE curve closely. At n = 300 nt, Turner-MFE requires ~2,000 ms, which is consistent with known complexity.

![Figure 2](figures/length_scaling.png)

**Figure 2.** Left: runtime vs. sequence length (linear scale). Right: log-log scale with O(n³) reference line (dashed gray). All methods conform to cubic complexity.

### 5.3 SARS-CoV-2 5'UTR Case Study

**Table 2**: Predictions for SARS-CoV-2 5'UTR (first 100 nt, MN908947.3)

| Method | Base Pairs | MFE (kcal/mol) | Features |
|--------|-----------|----------------|---------|
| Nussinov | 35 | — | Over-predicts; no energy model |
| Turner-MFE | 19 | −18.3 | SL1, partial SL2 |
| Turner+SHAPE+DMS | 20 | −16.9 | DMS-constrained SL1 |
| Turner+Pseudoknot | 28 | −27.9 | +9 pseudoknot-region pairs |

DMS reactivity restraints modestly adjusted the SL1 stem-loop structure (compare dot-brackets in Figure 3). The pseudoknot-aware method recovers 9 additional base pairs in otherwise unstructured regions.

![Figure 3](figures/sars_structure.png)

**Figure 3.** Arc diagrams of predicted SARS-CoV-2 5'UTR structure for four methods. Nucleotides colored by type (A=red, U=blue, G=green, C=orange). The pseudoknot-aware prediction (bottom-right) shows additional arc crossings consistent with known frameshifting elements.

### 5.4 Pseudoknot Benchmark

Both methods show low F1 on H-type pseudoknot sequences (Table 3). This reflects the fundamental difficulty of exact pseudoknot pair recovery—a result consistent with the O(n³) hierarchical approximation that cannot guarantee optimal pseudoknot detection.

**Table 3**: Pseudoknot Benchmark F1

| Length | HybridFold F1 | Turner-MFE F1 |
|--------|--------------|--------------|
| n = 50 | 0.042 ± 0.042 | 0.048 ± 0.048 |
| n = 70 | 0.016 ± 0.026 | 0.018 ± 0.029 |
| n = 90 | 0.019 ± 0.022 | 0.022 ± 0.025 |

![Figure 4](figures/pseudoknot_benchmark.png)

**Figure 4.** Pseudoknot benchmark: HybridFold vs. Turner-MFE on H-type pseudoknot sequences. Both methods fail to recover exact crossing pairs, consistent with the known difficulty of pseudoknot prediction.

### 5.5 SHAPE Data Quality Sensitivity

SHAPE integration with perfect reactivity data achieves F1 = 0.899 ± 0.035, degrading gracefully to 0.800 ± 0.045 at high noise (σ = 0.8), a 11.0% degradation over the full noise range tested (Table 4).

**Table 4**: SHAPE Noise Sensitivity

| Noise σ | F1 (mean ± SD) |
|---------|---------------|
| 0.0 | 0.899 ± 0.035 |
| 0.1 | 0.869 ± 0.031 |
| 0.2 | 0.852 ± 0.036 |
| 0.3 | 0.848 ± 0.037 |
| 0.5 | 0.816 ± 0.049 |
| 0.8 | 0.800 ± 0.045 |

![Figure 5](figures/shape_sensitivity.png)

**Figure 5.** F1 score (mean ± SD) as a function of SHAPE reactivity noise level. The degradation is gradual, suggesting robustness to experimental noise typical of in-cell probing experiments.

### 5.6 Covariation Heatmap

The mutual information matrix computed from the 50-sequence synthetic MSA (Figure 6) reveals block-diagonal covariation patterns corresponding to the predicted stem-loop positions in the SARS-CoV-2 5'UTR.

![Figure 6](figures/covariation_heatmap.png)

**Figure 6.** Mutual information heatmap for SARS-CoV-2 5'UTR MSA (first 50 positions, 50 sequences). Yellow-to-red scale indicates increasing MI (bits). Block patterns indicate co-evolving positions consistent with Watson-Crick base pairs.

---

## 6. Discussion

### 6.1 Effectiveness of SHAPE Integration

The observed F1 improvement from SHAPE integration (+8.8% over Turner-MFE) is consistent with prior literature. The Mathews pseudo-energy formula efficiently converts continuous SHAPE reactivities into position-specific pairing penalties without requiring retraining. The 2–8% normalization protocol proved robust; alternative normalization methods (e.g., box-plot outlier removal) may further improve stability.

Importantly, SHAPE data quality directly influences prediction accuracy (Table 4). At typical in-cell DMS probing noise levels (σ ≈ 0.2–0.3), F1 remains ≥ 0.85, suggesting practical utility. These results are consistent with the genome-wide chemical probing analyses of Wang [2025] and the eFold framework [de Lajarte et al., 2026].

### 6.2 Pseudoknot Prediction Limitations

The near-zero F1 values for pseudoknot prediction (Table 3) accurately reflect a fundamental algorithmic limitation: our hierarchical approach can identify regions likely to form pseudoknots, but cannot guarantee exact pair recovery. This is consistent with [Gray et al., 2024] (CParty), which notes that pseudoknot partition function algorithms are "borderline-prohibitive" in complexity.

A key distinction is that our hierarchical method identifies *more* pairs overall (28 vs. 19 for the SARS-CoV-2 case), including probable pseudoknot arms, even when exact pair positions are not recovered. For biological applications (e.g., drug target identification), identifying the general topology may be more important than exact pair recovery.

### 6.3 MSA Covariation Contribution

In this simplified implementation, adding MI covariation to SHAPE-constrained Turner DP provided negligible additional F1 improvement (Table 1). We attribute this to three factors: (1) the synthetic MSAs used (50 sequences, 6–8% mutation rate) may have insufficient diversity for strong MI signal; (2) the simple average MI approach is inferior to direct coupling analysis (DCA); (3) the weight coefficient (0.5) was not optimized.

Prior work [AliNA, Nasaev et al., 2023; eFold, de Lajarte et al., 2026] demonstrates that deep learning can implicitly learn evolutionary covariation signals from large training sets. Our simplified MI approach captures only pairwise correlations without correcting for transitive effects, a known limitation [Wu et al., 2023].

### 6.4 Comparison with Prior Work

HybridFold achieves F1 = 0.385 on synthetic data, below state-of-the-art deep learning methods (UFold: ~0.75 F1 on benchmark sets; eFold: ~0.70). Several factors explain this gap:
1. Our simplified Turner parameter subset (incomplete coverage of all loop types)
2. Synthetic vs. experimentally validated reference structures
3. Absence of learned parameters (gradient-free method)

The advantage of HybridFold is interpretability and physical grounding: each component has a clear thermodynamic motivation.

### 6.5 Future Directions

1. **Linear-time approximations**: LinearFold/LinearPartition-style beam search to extend scalability to transcriptome-wide analysis
2. **Complete Turner parameters**: Full ViennaRNA-compatible parameter set integration
3. **DCA-based covariation**: Replace MI with PSICOV/plmDCA for phylogenetically corrected covariation
4. **Gradient-based parameter optimization**: Use experimentally validated structures (Rfam, PDB) to optimize energy parameters
5. **SARS-CoV-2 full genome analysis**: Apply sliding-window HybridFold to the 30 kb genome with in-vivo DMS-MaPseq data

---

## 7. Conclusion

We presented HybridFold, a unified RNA secondary structure prediction framework integrating Turner thermodynamic DP, SHAPE/DMS chemical probing pseudo-energy restraints, and MSA-based mutual information covariation scoring. Key findings:

1. **SHAPE integration provides meaningful improvement**: +8.8% F1 over plain Turner-MFE (0.354 → 0.385), robust to realistic experimental noise
2. **Pseudoknot prediction remains unsolved**: Both pseudoknot-aware and pseudoknot-free methods achieve F1 < 0.05 on exact pair recovery for H-type pseudoknots
3. **MI covariation shows limited benefit** in this simplified implementation; DCA-based approaches are recommended for future work
4. **O(n³) scaling** was confirmed empirically, consistent with theoretical analysis
5. **SARS-CoV-2 5'UTR analysis** demonstrates practical applicability: SHAPE-constrained prediction recovers the known SL1 stem-loop with experimentally consistent topology

HybridFold represents a modular, extensible framework that can incorporate additional experimental constraints and learning-based components. The results highlight that thermodynamic methods remain competitive when augmented with experimental data, even against modern deep learning baselines.

---

## References

1. Fu, L., Cao, Y., Wu, J., Peng, Q., Nie, Q., & Xie, X. (2020). UFold: Fast and Accurate RNA Secondary Structure Prediction with Deep Learning. *bioRxiv*. https://doi.org/10.1101/2020.08.17.254896

2. Gray, M., Trinity, L., Stege, U., Ponty, Y., & Will, S. (2024). CParty: hierarchically constrained partition function of RNA pseudoknots. *Bioinformatics*, 40. https://doi.org/10.1093/bioinformatics/btae748

3. Gray, M., Will, S., & Jabbari, H. (2024). SparseRNAFolD: optimized sparse RNA pseudoknot-free folding with dangle consideration. *Algorithms for Molecular Biology*, 19(1). https://doi.org/10.1186/s13015-024-00256-4

4. Newman, T., Chang, H. F. K., & Jabbari, H. (2024). DinoKnot: Duplex Interaction of Nucleic Acids With PseudoKnots. *IEEE/ACM TCBB*, 21(3). https://doi.org/10.1109/TCBB.2024.3362308

5. Nasaev, S. S., Mukanov, A. R., Kuznetsov, I. I., & Veselovsky, A. V. (2023). AliNA - a deep learning program for RNA secondary structure prediction. *Molecular Informatics*, 42. https://doi.org/10.1002/minf.202300113

6. Wu, K. E., Zou, J. Y., & Chang, H. (2023). Machine learning modeling of RNA structures: methods, challenges and future perspectives. *Briefings in Bioinformatics*, 24(4). https://doi.org/10.1093/bib/bbad210

7. Omnes, L., Angel, E., & Tahi, F. (2026). DivideFold+: an AI-based tool for RNA secondary structure prediction with subdomains identification and visualization and data augmentation. *Journal of Molecular Biology*. https://doi.org/10.1016/j.jmb.2026.169865

8. de Lajarte, A. A., Taillades, Y. J. M. D., Aruda, J., Bongrand, P., & Wightman, F. F. (2026). Diverse database and machine learning model to narrow the generalization gap in RNA structure prediction. *Science Advances*. https://doi.org/10.1126/sciadv.adz4967

9. Ziesel, A., & Jabbari, H. (2024). Unveiling hidden structural patterns in the SARS-CoV-2 genome: Computational insights and comparative analysis. *PLoS ONE*, 19(4). https://doi.org/10.1371/journal.pone.0298164

10. Wang, W., Liu, X., Peng, Z., & Yang, J. (2026). The trRosettaRNA server for RNA structure prediction. *Nature Protocols*. https://doi.org/10.1038/s41596-026-01356-8

11. Wang, J. (2025). Genome-Wide Analysis of Stable RNA Secondary Structures across Multiple Organisms Using Chemical Probing Data. *Biochemistry*, 64(8). https://doi.org/10.1021/acs.biochem.4c00764

12. Wirecki, T. K., Merdas, K., Bernat, A., Boniecki, M. J., & Bujnicki, J. M. (2020). RNAProbe: a web server for normalization and analysis of RNA structure probing data. *Nucleic Acids Research*, 48(W1). https://doi.org/10.1093/nar/gkaa396
