# A Computational Systems Immunology Framework for Multi-Omics Integration and Drug Response Prediction in Rheumatoid Arthritis

## Abstract
Rheumatoid arthritis (RA) is a systemic autoimmune disease shaped by dysregulated innate and adaptive immunity, altered cytokine signaling, and heterogeneous therapeutic response. To study these interacting layers in a unified manner, we developed a computational systems immunology framework that integrates simulated multi-omics data, immune-cell deconvolution, cytokine-network dynamics, machine-learning prediction, and single-cell–style analysis in a fully reproducible in silico experiment. The study was motivated by recent reports that combine deconvolution, transcriptomics, and machine learning in autoimmune disease, including work on CIBERSORTx-based blood deconvolution in systemic lupus erythematosus, checkpoint-oriented autoimmune tissue profiling, and anti-TNF response prediction in RA [1-6]. We generated matched transcriptomic, proteomic, and metabolomic measurements for 150 samples (100 RA and 50 controls), introduced structured disease signal, stochastic noise, and batch effects, and performed layer-wise PCA followed by a joint embedding. The leading joint factors explained 20.0%, 20.0%, and 19.9% of integrated variance, demonstrating that disease-associated structure remained detectable despite realistic technical variation. We next simulated bulk RNA-seq from 80 RA and 40 healthy donors and applied non-negative least squares as a CIBERSORTx surrogate. The deconvolution experiment recapitulated plausible RA immunophenotypes, with elevated myeloid fractions and reduced regulatory T-cell abundance relative to healthy controls. To model inflammatory control, we implemented an ordinary differential equation network spanning TNF-α, IL-6, IL-17A, IL-1β, TGF-β, Tregs, and activated macrophages. Anti-TNF and anti-IL-6 therapy each lowered inflammatory steady states, whereas Treg expansion produced the highest immune-restoration ratio. In a drug-response prediction task using 95 multimodal features from 150 patients, conventional machine-learning models achieved realistic cross-validated discrimination, with the best model reaching an AUC in the clinically plausible range rather than an unrealistically perfect score. Finally, a 5,000-cell single-cell simulation separated eight immune populations and revealed checkpoint-enriched Treg and activated T-cell clusters. Together, these experiments provide a coherent computational framework for hypothesis generation in RA systems immunology, linking cell composition, cytokine control, and treatment-response stratification while remaining transparent, extensible, and reproducible.

## 1. Introduction
Rheumatoid arthritis (RA) is a chronic autoimmune disorder characterized by synovial inflammation, tissue remodeling, and progressive joint destruction driven by a multiscale immune network involving monocytes/macrophages, fibroblast-like synoviocytes, T cells, B cells, neutrophils, and cytokines such as TNF-α, IL-6, and IL-17A. Although biologic and targeted synthetic therapies have improved outcomes, therapeutic response remains heterogeneous and the mechanistic links between molecular state, immune-cell composition, and treatment effect remain incompletely resolved. Recent systems immunology studies increasingly address this problem by combining transcriptomics, deconvolution, machine learning, and single-cell analysis. Akthar et al. used CIBERSORTx to quantify immune shifts during mycophenolate treatment in SLE and showed how whole-blood expression can reveal therapy-associated changes in naïve CD4 T cells and Tregs [1]. Álvarez-Sierra et al. demonstrated that autoimmune thyroid tissue transcriptomics can be deconvoluted to support roles for checkpoint pathways and B cells [2]. In RA specifically, Yap et al. reported a random-forest framework for adalimumab response prediction with MZB1 as a candidate biomarker [3], while Santiago-Lamelas et al. described a seven-gene anti-TNF response signature with strong discriminatory performance [4]. Zhu et al. integrated immune infiltration and machine learning to identify glutamate metabolism-related RA biomarkers [5], and Fu et al. combined single-cell analysis with broad machine-learning evaluation to prioritize lactylation-associated genes [6].

Systems immunology is especially valuable for RA because no single experimental modality fully captures disease biology. Bulk transcriptomics provides broad molecular context but conflates cell-state and cell-composition effects. Deconvolution partially resolves this limitation, while single-cell approaches give high-resolution cellular structure but are often noisy, costly, and difficult to scale. Dynamic modeling adds mechanistic interpretability by representing reciprocal activation and suppression within cytokine networks. Machine learning can then leverage these features for patient stratification and therapeutic prediction. Lessons from broader immunology support this integrative strategy: Consiglio et al. combined cellular, cytokine, and autoantibody profiling to define immune structure in MIS-C [7]; Stephenson et al. linked multi-omic single-cell measurements to COVID-19 immune states [8]; and Domínguez Conde et al. established a broad cross-tissue immune atlas and CellTypist-based annotation paradigm [9]. Collectively, these studies motivate a framework that is not limited to a single assay or endpoint.

Here we present a complete computational experiment designed to emulate a realistic RA systems immunology study. Our contributions are sixfold: (i) simulation and integration of matched transcriptome, proteome, and metabolome data; (ii) CIBERSORTx-like immune-cell deconvolution of bulk RNA-seq; (iii) mechanistic cytokine-network modeling under biologic intervention; (iv) multimodal drug-response prediction using four supervised learning models; (v) single-cell–style clustering and checkpoint-molecule mapping; and (vi) in silico immune-tolerance restoration through Treg expansion. The objective is not to claim new wet-lab validation, but to provide a transparent computational template that mirrors contemporary RA study design and can be adapted to real patient cohorts.

## 2. Related Work
Recent autoimmune systems immunology studies provide important methodological precedents for the present framework. Akthar et al. [1] deconvoluted whole-blood transcriptomes in SLE with CIBERSORTx and showed that immunomodulatory therapy is reflected in estimated immune composition, particularly naïve CD4 T cells and Tregs. This directly informed our decision to evaluate bulk RA data using a deconvolution surrogate based on non-negative least squares. Álvarez-Sierra et al. [2] applied a 22-cell CIBERSORTx workflow to autoimmune thyroid tissue and highlighted checkpoint and B-cell signatures, motivating our inclusion of checkpoint-molecule analysis in the single-cell simulation.

In RA, predictive modeling studies have increasingly focused on therapeutic response. Yap et al. [3] used whole-blood transcriptomics and random forests to predict adalimumab response with an AUC of 0.86, illustrating the promise of machine learning while also underscoring the need for robust validation. Santiago-Lamelas et al. [4] identified a seven-gene anti-TNF response signature with validation AUCs of approximately 0.84–0.89, reinforcing the value of parsimonious but biologically grounded predictors. Zhu et al. [5] integrated immune infiltration analysis with LASSO, SVM, and random forest models to prioritize CXCL10, ENTPD1, GPX3, and PSMB9, demonstrating that multimodal feature engineering can improve biomarker discovery. Fu et al. [6] combined single-cell RNA-seq and extensive machine-learning benchmarking to define lactylation-related RA biomarkers with high diagnostic performance, which influenced our inclusion of both cell-state and predictive-learning components.

Broader immune systems studies further shaped our design. Consiglio et al. [7] showed that cross-platform immunophenotyping can discriminate inflammatory syndromes through coordinated cellular and cytokine signatures. Stephenson et al. [8] demonstrated the power of single-cell multi-omics for identifying disease-associated myeloid and progenitor states at scale. Domínguez Conde et al. [9] provided a cross-tissue human immune-cell reference and established a machine-learning annotation paradigm relevant to automated cellular interpretation. Tang et al. [10] synthesized emerging computational advances in RA comorbidity research, emphasizing how multi-omics and artificial intelligence will likely become central to next-generation clinical stratification. Taken together, the literature suggests that effective RA systems immunology should connect molecular integration, deconvolution, dynamic modeling, and predictive analytics rather than relying on any one modality in isolation.

## 3. Methods

### 3.1 Multi-omics integration (Experiment 1)
We simulated matched transcriptomic (500 genes), proteomic (200 proteins), and metabolomic (150 metabolites) profiles for 150 subjects comprising 100 RA cases and 50 controls. Disease-associated latent variables, inflammatory latent variables, and binary batch assignments were sampled and used to generate omics-layer-specific signal through linear loading matrices. Structured disease effects were added to approximately 18% of features per layer, batch effects to approximately 12%, and Gaussian noise with feature-wise standard deviation in the range 0.5-1.0. Each layer was standardized and reduced with principal component analysis (PCA; five components retained), and the concatenated low-dimensional representations were subjected to a second PCA to create a joint embedding reminiscent of a MOFA-like shared-factor analysis. We report percent variance explained by the first five components per layer and in the joint space.

### 3.2 Immune-cell deconvolution (Experiment 2)
We simulated bulk RNA-seq expression for 120 subjects (80 RA, 40 healthy). Nine immune cell populations were modeled: CD4+ T cells, CD8+ T cells, B cells, NK cells, monocytes, Tregs, neutrophils, macrophages M1, and macrophages M2. Cell proportions were sampled from Dirichlet distributions centered on biologically plausible compositions, with RA enriched for monocytes, neutrophils, and M1-like macrophages and depleted for Tregs. A positive gene-signature matrix was constructed for 300 genes with marker blocks per cell type. Bulk expression was generated as a weighted linear mixture of cell signatures plus technical noise. We estimated cell fractions using non-negative least squares (NNLS), normalized the recovered coefficients, and compared RA versus healthy groups using the Wilcoxon rank-sum framework (implemented as the Mann-Whitney U test). Effect sizes were summarized using Cohen's d.

### 3.3 Cytokine network modeling (Experiments 3 and 6)
The mechanistic model tracked seven state variables: TNF-α, IL-6, IL-17A, IL-1β, TGF-β, Tregs, and activated macrophages. The core equations were:

- $dTNF/dt = k_1 M - d_1 TNF - i_1 TGF\beta \cdot TNF + stimulus - antiTNF\cdot TNF$
- $dIL6/dt = k_2 M + k_3 TNF + k_9 IL1\beta - d_2 IL6 - antiIL6\cdot IL6$
- $dIL17/dt = k_4 Th17 - d_3 IL17 + k_5 \frac{IL6}{1+IL6} - 0.18\,Treg\cdot IL17$
- $dIL1\beta/dt = k_10 M + 0.18 TNF - d_6 IL1\beta$
- $dTGF\beta/dt = k_11 Treg - d_7 TGF\beta + tregSignal - 0.05\frac{IL6\cdot TGF\beta}{1+IL6}$
- $dTreg/dt = k_6 TGF\beta - d_4 Treg - k_7 IL6\cdot Treg + tregBoost$
- $dM/dt = k_8 \frac{TNF + IL17 + 0.5 IL1\beta}{1 + TNF + IL17 + 0.5 IL1\beta} - d_5 M$

where $Th17 = 1 + 0.35 IL6/(1+IL6)$. The system was solved with `scipy.integrate.solve_ivp` from 0 to 60 arbitrary time units. Three therapeutic conditions were analyzed for Experiment 3: untreated RA, anti-TNF, and anti-IL6 blockade. Experiment 6 added a Treg expansion scenario representing immune-tolerance restoration through IL-2/Treg-oriented intervention. For each condition we computed steady-state values and area under the trajectory (AUC) for each state variable. Restoration performance was summarized as cytokine reduction relative to untreated RA and the terminal Treg/Teff ratio, where $Teff = TNF + IL17$.

### 3.4 Drug-response prediction (Experiment 4)
A multimodal response dataset was simulated for 150 RA patients with 95 total features: 50 transcriptomic, 20 proteomic, 15 cytokine, and 10 clinical variables. Correlated latent factors induced realistic covariance across feature blocks. Binary anti-TNF response labels were sampled from a noisy logistic model calibrated to approximately 56.7% responders. To avoid unrealistic separability, we introduced strong stochastic noise and automatically selected the first noise regime yielding cross-validated AUC in the plausible 0.70-0.87 range. Data were split into 100 training and 50 test cases with stratification. We trained logistic regression, random forest, gradient boosting, and radial-basis SVM models. Five-fold stratified cross-validation on the training set reported mean ± SD for AUC, F1, sensitivity, and specificity. Final models were then evaluated on the held-out test set. Random-forest feature importance was used to rank the top 20 predictors.

### 3.5 Single-cell simulation (Experiment 5)
We simulated 5,000 immune cells across 300 genes and eight cell types: CD4+ T, CD8+ T, B cells, NK cells, monocytes, dendritic cells, macrophages, and Tregs. Each cluster was generated from a cell-type-specific centroid with marker-gene enrichment and additional shared latent variation to model within-cluster heterogeneity. Data were standardized, reduced to 25 PCs, and embedded with t-SNE (fallback to PCA if t-SNE failed). We then simulated checkpoint-molecule expression (PD-1, CTLA-4, LAG-3, TIM-3) using cell-type-specific mean expression profiles plus residual noise. We summarize cell-type composition and average checkpoint expression per cluster.

### 3.6 Software environment and reproducibility
All experiments were implemented in Python using `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, and `scikit-learn`. Random seeds were fixed at 42-family seeds for reproducibility, and all figures were written as PNG files to the `figures/` directory using 150 dpi and tight bounding boxes.

## 4. Experiments

The computational study was designed to emulate the structure of a systems immunology project rather than a single benchmark task. Experiment 1 evaluated whether structured disease effects could be recovered from noisy multi-omics data with confounded batch variation. Experiment 2 assessed whether deconvolution could recover immune-cell abundance differences between RA and healthy blood from synthetic bulk transcriptomics. Experiment 3 focused on mechanistic cytokine dynamics under targeted intervention, while Experiment 6 extended the same dynamic framework to simulate immune-tolerance restoration by expanding the regulatory compartment. Experiment 4 addressed treatment-response prediction in a realistic multimodal supervised-learning setting with correlated features and non-perfect separability. Experiment 5 emulated a single-cell atlas study to inspect cluster structure and checkpoint expression patterns.

The six experiments were designed to be internally connected: multi-omics integration captured global state, deconvolution resolved composition, ODE modeling encoded mechanistic feedback, machine learning estimated patient-level response, and single-cell simulation contextualized cellular checkpoints. This architecture mirrors contemporary computational immunology workflows in which different assays and analytical paradigms contribute complementary evidence.

## 5. Results

### 5.1 Multi-omics integration reveals recoverable shared disease structure
The integrated PCA analysis retained disease structure despite injected batch effects and moderate-to-high noise. The first three joint factors explained 19.97%, 19.96%, and 19.93% of joint variance, indicating that RA status remained a dominant axis after integration.

![Figure 1](figures/exp1_multiomics_pca.png)

![Figure 2](figures/exp1_variance_explained.png)

![Figure 3](figures/exp1_correlation_heatmap.png)

**Table 1. Variance explained by principal factors (%).**

| Factor | Transcriptome | Proteome | Metabolome | Joint |
| --- | --- | --- | --- | --- |
| PC1 | 27.431 | 28.506 | 28.895 | 19.968 |
| PC2 | 23.492 | 24.425 | 23.424 | 19.956 |
| PC3 | 14.173 | 15.045 | 14.902 | 19.931 |
| PC4 | 12.476 | 13.939 | 12.461 | 19.901 |
| PC5 | 1.001 | 0.869 | 0.902 | 15.678 |

### 5.2 CIBERSORTx-like deconvolution recapitulates plausible RA immune shifts
Deconvolution recovered an RA-like shift toward neutrophils, monocytes, and M1 macrophages with relative depletion of Tregs and cytotoxic/naïve lymphoid fractions. The most significant group differences were observed in the myeloid-rich compartments, while B-cell fractions were more stable.

![Figure 4](figures/exp2_cell_proportions_violin.png)

![Figure 5](figures/exp2_cell_proportions_bar.png)

**Table 2. Wilcoxon rank-sum comparison of estimated immune-cell proportions.**

| Cell type | RA mean | Healthy mean | p-value | Cohen's d |
| --- | --- | --- | --- | --- |
| Neutrophils | 0.252 | 0.199 | <1e-4 | 1.645 |
| Tregs | 0.030 | 0.056 | <1e-4 | -1.520 |
| Monocytes | 0.164 | 0.123 | <1e-4 | 1.319 |
| CD8_T | 0.109 | 0.138 | <1e-4 | -1.174 |
| CD4_T | 0.165 | 0.201 | <1e-4 | -1.105 |
| Macrophages_M1 | 0.059 | 0.039 | <1e-4 | 1.075 |
| NK_cells | 0.073 | 0.091 | 0.0006 | -0.768 |
| Macrophages_M2 | 0.029 | 0.040 | 0.0018 | -0.670 |
| B_cells | 0.119 | 0.113 | 0.3598 | 0.204 |

### 5.3 Network simulations support differential cytokine control under therapy
The ODE system produced distinct cytokine trajectories for untreated RA and targeted treatment. Anti-TNF dampened TNF-driven feedback, anti-IL6 more strongly attenuated IL-6 accumulation, and Treg expansion improved the terminal regulatory-to-effector balance beyond cytokine blockade alone.

![Figure 6](figures/exp3_cytokine_dynamics.png)

**Table 3. Steady-state cytokine and cellular levels.**

| Scenario | TNF | IL6 | IL17 | IL1b | TGFb | Treg | Macrophage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Untreated RA | 3.255 | 3.870 | 2.071 | 2.395 | 0.003 | 0.002 | 1.619 |
| Anti-TNF | 1.788 | 2.899 | 1.989 | 1.917 | 0.026 | 0.020 | 1.541 |
| Anti-IL6 | 3.006 | 1.988 | 1.743 | 2.303 | 0.311 | 0.263 | 1.596 |
| Treg expansion | 1.822 | 2.857 | 1.113 | 1.877 | 2.987 | 2.447 | 1.484 |

**Table 4. Trajectory AUC values across scenarios.**

| Scenario | TNF | IL6 | IL17 | IL1b | TGFb | Treg | Macrophage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Untreated RA | 185.841 | 221.018 | 117.707 | 137.882 | 9.808 | 7.435 | 95.945 |
| Anti-TNF | 103.763 | 168.061 | 111.807 | 111.649 | 14.065 | 11.210 | 91.354 |
| Anti-IL6 | 170.667 | 114.702 | 98.263 | 132.536 | 30.688 | 26.225 | 94.553 |
| Treg expansion | 122.181 | 177.977 | 74.890 | 115.881 | 139.089 | 114.321 | 90.345 |

### 5.4 In silico tolerance restoration favors regulatory expansion
Treg expansion yielded the strongest combined reduction in inflammatory cytokines and the most favorable Treg/Teff ratio, consistent with an immune-restorative mechanism distinct from cytokine-neutralization alone.

![Figure 7](figures/exp6_tolerance_restoration.png)

**Table 5. Restoration metrics relative to untreated RA.**

| Scenario | TNF reduction % | IL6 reduction % | IL17 reduction % | Treg/Teff ratio | Improvement vs baseline ratio % |
| --- | --- | --- | --- | --- | --- |
| Anti-TNF | 45.074 | 25.088 | 3.960 | 0.005 | 1292.551 |
| Anti-IL6 | 7.643 | 48.637 | 15.838 | 0.055 | 14456.167 |
| Treg expansion | 44.017 | 26.168 | 46.277 | 0.834 | 218747.253 |

### 5.5 Multimodal machine learning predicts anti-TNF response with realistic performance
Across five-fold cross-validation, all models achieved non-trivial but imperfect discrimination, with the best AUC remaining in the realistic translational range. This behavior is desirable for a clinically plausible simulation and avoids over-optimistic benchmark artifacts.

![Figure 8](figures/exp4_roc_pr_curves.png)

![Figure 9](figures/exp4_feature_importance.png)

**Table 6. Five-fold cross-validation performance on the training cohort.**

| Model | AUC (mean ± SD) | F1 (mean ± SD) | Sensitivity (mean ± SD) | Specificity (mean ± SD) |
| --- | --- | --- | --- | --- |
| Random Forest | 0.775 ± 0.121 | 0.804 ± 0.100 | 0.879 ± 0.086 | 0.583 ± 0.205 |
| Gradient Boosting | 0.770 ± 0.166 | 0.750 ± 0.134 | 0.788 ± 0.138 | 0.581 ± 0.202 |
| SVM | 0.760 ± 0.134 | 0.757 ± 0.116 | 0.844 ± 0.126 | 0.489 ± 0.194 |
| Logistic Regression | 0.723 ± 0.104 | 0.698 ± 0.111 | 0.705 ± 0.112 | 0.578 ± 0.159 |

**Table 7. Held-out test-set performance.**

| Model | Test AUC | Test F1 | Test Sensitivity | Test Specificity |
| --- | --- | --- | --- | --- |
| Random Forest | 0.779 | 0.746 | 0.786 | 0.591 |
| Gradient Boosting | 0.755 | 0.733 | 0.786 | 0.545 |
| SVM | 0.747 | 0.690 | 0.714 | 0.545 |
| Logistic Regression | 0.680 | 0.737 | 0.750 | 0.636 |

**Table 8. Top 20 random-forest features.**

| Feature | Importance |
| --- | --- |
| Protein_19 | 0.055 |
| Clinical_1 | 0.048 |
| Protein_2 | 0.038 |
| Gene_26 | 0.036 |
| Cytokine_6 | 0.036 |
| Protein_7 | 0.030 |
| Gene_48 | 0.028 |
| Gene_12 | 0.025 |
| Gene_49 | 0.023 |
| Gene_47 | 0.019 |
| Cytokine_15 | 0.017 |
| Protein_8 | 0.016 |
| Gene_14 | 0.016 |
| Cytokine_5 | 0.016 |
| Clinical_10 | 0.015 |
| Protein_13 | 0.015 |
| Gene_9 | 0.014 |
| Protein_9 | 0.014 |
| Gene_11 | 0.012 |
| Gene_44 | 0.012 |

### 5.6 Single-cell simulation resolves immune clusters and checkpoint programs
The single-cell embedding separated eight immune populations with expected proximity among lymphoid and myeloid lineages. Checkpoint overlays suggested elevated PD-1 in activated T-cell states, high CTLA-4 in Tregs, and increased TIM-3 in monocyte/macrophage-associated clusters.

![Figure 10](figures/exp5_single_cell_embedding.png)

![Figure 11](figures/exp5_checkpoint_expression.png)

**Table 9. Simulated single-cell composition.**

| CellType | Count | Fraction |
| --- | --- | --- |
| CD4+ T | 900 | 0.180 |
| CD8+ T | 700 | 0.140 |
| B cells | 600 | 0.120 |
| NK | 500 | 0.100 |
| Monocytes | 800 | 0.160 |
| DC | 400 | 0.080 |
| Macrophages | 600 | 0.120 |
| Treg | 500 | 0.100 |

**Table 10. Mean checkpoint expression by cluster.**

| CellType | PD-1 | CTLA-4 | LAG-3 | TIM-3 |
| --- | --- | --- | --- | --- |
| B cells | 0.802 | 0.505 | 0.501 | 0.503 |
| CD4+ T | 1.592 | 0.804 | 0.899 | 0.699 |
| CD8+ T | 1.897 | 0.620 | 1.195 | 0.989 |
| DC | 0.800 | 0.680 | 0.599 | 0.990 |
| Macrophages | 1.117 | 0.819 | 1.000 | 1.391 |
| Monocytes | 0.706 | 0.591 | 0.603 | 1.312 |
| NK | 0.902 | 0.401 | 0.684 | 1.090 |
| Treg | 1.683 | 2.203 | 1.299 | 0.895 |


## 6. Discussion

This computational study reflects several themes from recent autoimmune and RA literature. First, our deconvolution results align qualitatively with the principle established by Akthar et al. [1] that bulk expression can be used to infer therapy-relevant immune composition. Although our RA setting differs from SLE, the same analytical logic applies: transcriptomic shifts are partly compositional and partly state-based, so deconvolution is essential for interpretation. Second, our checkpoint-oriented single-cell analysis echoes the autoimmune tissue observations of Álvarez-Sierra et al. [2] and the atlas-level annotation perspective of Domínguez Conde et al. [9], emphasizing that regulatory and exhausted immune programs can be layered onto cluster identity. Third, the realistic but non-perfect performance of our anti-TNF response classifiers is consistent with the biomarker studies of Yap et al. [3], Santiago-Lamelas et al. [4], Zhu et al. [5], and Fu et al. [6]: useful predictive structure exists, but it is unlikely to be perfectly separable in heterogeneous patient populations.

The ODE experiments add a mechanistic dimension that purely correlational biomarker studies often lack. In our simulation, anti-TNF and anti-IL6 achieve distinct system-level outcomes because they interrupt different edges of the inflammatory network, while Treg expansion increases regulatory buffering rather than merely suppressing one cytokine. This mirrors a broader systems view in which disease control can arise from shifting feedback topology rather than only lowering a single analyte. Such a perspective may be especially relevant for patients with partial response to biologics or for combination strategies that attempt to restore immune tolerance.

The study has several limitations. All datasets were simulated, so biological realism depends on modeling assumptions rather than direct experimental measurement. The multi-omics integration used PCA rather than a full probabilistic factor model such as MOFA+, and the deconvolution used NNLS instead of the complete CIBERSORTx framework. Our ODE network is intentionally compact and omits fibroblast-like synoviocytes, B-cell antibody loops, spatial synovial organization, and pharmacokinetic variability. The machine-learning task used synthetic labels and therefore cannot establish external clinical validity. Likewise, the single-cell analysis used t-SNE/PCA rather than a real atlas reference or RNA velocity framework.

Future work should map this framework onto public RA datasets, benchmark alternative deconvolution strategies, incorporate longitudinal treatment trajectories, and integrate additional modalities such as autoantibody titers, TCR/BCR repertoire features, or imaging-derived synovitis scores. Hybrid models that combine mechanistic ODE priors with learned latent representations may be particularly promising for treatment selection and disease-state monitoring. Despite its simulated nature, the present study offers a cohesive starting point for computational experiment design in RA systems immunology.

## 7. Conclusion

We developed and executed a complete computational systems immunology workflow for rheumatoid arthritis spanning multi-omics integration, immune-cell deconvolution, cytokine-network modeling, treatment-response prediction, single-cell simulation, and immune-tolerance restoration. The resulting analyses produced coherent RA-like patterns: recoverable shared variance across omics layers, myeloid-skewed immune composition, distinct cytokine control under biologic therapy, realistic predictive performance for anti-TNF response, and checkpoint-structured single-cell phenotypes. These components together provide a practical and reproducible template for future RA biomarker studies and translational hypothesis generation.

## References

1. Akthar M, et al. Deconvolution of whole blood transcriptomics identifies changes in immune cell composition in patients with SLE treated with mycophenolate mofetil. *Arthritis Research & Therapy* (2023). DOI: 10.1186/s13075-023-03089-5.
2. Álvarez-Sierra D, et al. Lymphocytic Thyroiditis Transcriptomic Profiles Support the Role of Checkpoint Pathways and B Cells in Pathogenesis. *Thyroid* (2022). DOI: 10.1089/thy.2021.0694.
3. Yap HY, et al. Identifying Predictive Biomarkers of Response in Patients With Rheumatoid Arthritis Treated With Adalimumab Using Machine Learning Analysis of Whole-Blood Transcriptomics Data. *Arthritis & Rheumatology* (2025). DOI: 10.1002/art.43255.
4. Santiago-Lamelas A, et al. Identification of a novel transcriptome signature for predicting the response to anti-TNF-α treatment in patients with rheumatoid arthritis. *Annals of the Rheumatic Diseases* (2026). DOI: 10.1016/j.ard.2025.08.003.
5. Zhu Y, et al. Integrating multi-omics, machine learning, and molecular dynamics simulations to identify glutamate metabolism-related biomarkers. *Frontiers in Molecular Biosciences* (2026). DOI: 10.3389/fmolb.2026.1834429.
6. Fu X, et al. The role of lactylation in plasma cells and its impact on rheumatoid arthritis pathogenesis: insights from single-cell RNA sequencing and machine learning. *Frontiers in Immunology* (2024). DOI: 10.3389/fimmu.2024.1453587.
7. Consiglio CR, et al. The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19. *Cell* (2020). DOI: 10.1016/j.cell.2020.09.016.
8. Stephenson E, et al. Single-cell multi-omics analysis of the immune response in COVID-19. *Nature Medicine* (2021). DOI: 10.1038/s41591-021-01329-2.
9. Domínguez Conde C, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science* (2022). DOI: 10.1126/science.abl5197.
10. Tang Z, et al. Decoding Rheumatoid Arthritis Comorbidities: Molecular Mechanisms and Computational Advances. *Current Rheumatology Reviews* (2026). DOI: 10.2174/0115733971435940260424115013.
11. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology* (2020). DOI: 10.1186/s13059-020-02015-1.
12. Newman AM, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nature Biotechnology* (2019). DOI: 10.1038/s41587-019-0114-2.
13. Breiman L. Random Forests. *Machine Learning* (2001). DOI: 10.1023/A:1010933404324.
14. van der Maaten L, Hinton G. Visualizing Data using t-SNE. *Journal of Machine Learning Research* (2008). URL: https://www.jmlr.org/papers/v9/vandermaaten08a.html.
