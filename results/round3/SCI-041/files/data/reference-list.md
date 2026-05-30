# Reference List

[1] Rives A, Meier J, Sercu T, et al. (2021). **Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences**. *PNAS*. DOI: 10.1073/pnas.2016239118

   → ESM-1b; unsupervised protein LM; zero-shot mutation scoring

[2] Meier J, Rao R, Verkuil R, et al. (2021). **Language models enable zero-shot prediction of the effects of mutations on protein function**. *NeurIPS 2021 (bioRxiv)*. DOI: 10.1101/2021.07.09.450648

   → ESM-1v masked marginal scoring; SOTA zero-shot fitness

[3] Brandes N, Ofer D, Peleg Y, et al. (2022). **ProteinBERT: a universal deep-learning model of protein sequence and function**. *Bioinformatics*. DOI: 10.1093/bioinformatics/btac020

   → BERT + GO annotation; near-SOTA on 9 protein property benchmarks

[4] Lin Z, Akin H, Rao R, et al. (2022). **Evolutionary-scale prediction of atomic level protein structure with a language model**. *Science (bioRxiv)*. DOI: 10.1101/2022.07.20.500902

   → ESM-2 8M–15B; ESMFold; 60× faster than AlphaFold2

[5] Zeng S, Wang D, Jiang L, Xu D. (2024). **Parameter-efficient fine-tuning on large protein language models improves signal peptide prediction**. *Genome Research*. DOI: 10.1101/gr.279132.124

   → LoRA+ESM-2; 87.3% MCC gain in low-data regime vs full fine-tuning

[6] Notin P, Kollasch AW, Ritter DP, et al. (2023). **ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction**. *bioRxiv*. DOI: 10.1101/2023.12.07.570727

   → 250+ DMS assays; evaluation of 70+ models; zero-shot & supervised

[7] Gelman S, Johnson B, Freschlin CR, et al. (2025). **Biophysics-based protein language models for protein engineering**. *Nature Methods*. DOI: 10.1038/s41592-025-02776-2

   → METL: biophysical simulation + PLM fine-tuning; 64-example GFP design

[8] Ding K, Chin MA, Zhao Y, et al. (2024). **Machine learning-guided co-optimization of fitness and diversity in enzyme engineering**. *Nature Communications*. DOI: 10.1038/s41467-024-50698-y

   → MODIFY; ESM zero-shot + diversity; cytochrome c engineering

[9] Zhang Q, Chen W, Qin M, et al. (2025). **Integrating protein language models and automatic biofoundry for enhanced protein evolution**. *Nature Communications*. DOI: 10.1038/s41467-025-56751-8

   → ESM-2 + closed-loop biofoundry; 2.4-fold tRNA synthetase in 10 days

[10] Ferruz N, Schmidt S, Höcker B. (2022). **ProtGPT2 is a deep unsupervised language model for protein design**. *Nature Communications*. DOI: 10.1038/s41467-022-32007-7

   → GPT-2 on UniRef50; de novo well-folded protein generation

[11] Alley EC, Khimulya G, Biswas S, et al. (2019). **Unified rational protein engineering with sequence-based deep representation learning**. *Nature Methods*. DOI: 10.1038/s41592-019-0598-1

   → UniRep LSTM; unified representation for diverse protein engineering tasks

[12] Buehler EL, Buehler MJ. (2024). **X-LoRA: Mixture of low-rank adapter experts for large language models with applications in protein mechanics**. *APL Machine Learning*. DOI: 10.1063/5.0203126

   → MoE-LoRA dynamic gating; protein mechanics and molecular design

[13] Rao R, Liu J, Verkuil R, et al. (2021). **MSA Transformer**. *ICML 2021*. DOI: 10.1101/2021.02.12.430858

   → Row-column attention over MSA; unsupervised structure learning SOTA

