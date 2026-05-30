# Reference List — Epigenetic Clock Study

## Databases Searched

| Database | Tool | Status |
|----------|------|--------|
| PubMed | `PubMed_search_articles` | ✅ Success |
| PMC | (secondary) | ✅ Success |
| OpenAlex | `openalex_literature_search` | ⚠️ Returned off-topic results for keyword query — not used |
| ArXiv | `ArXiv_search_papers` | ❌ Network timeout — not included |

## Included References

1. (Horvath, 2013) Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115. PMID: 24138928. https://doi.org/10.1186/gb-2013-14-10-r115

2. (Lu, 2019) Lu, A. T., Quach, A., Wilson, J. G., Reiner, A. P., Aviv, A., et al. (2019). DNA methylation GrimAge strongly predicts lifespan and healthspan. *Aging*, 11(2), 303–327. PMID: 30669119. https://doi.org/10.18632/aging.101684

3. (Prosz, 2024) Prosz, A., Pipek, O., Börcsök, J., Palla, G., & Szallasi, Z. (2024). Biologically informed deep learning for explainable epigenetic clocks. *Scientific Reports*, 14, 1439. PMID: 38225268. https://doi.org/10.1038/s41598-023-50495-5

4. (Kalyakulina, 2025) Kalyakulina, A., Yusipov, I., Trukhanov, A., Franceschi, C., & Moskalev, A. (2025). EpInflammAge: Epigenetic-Inflammatory Clock for Disease-Associated Biological Aging Based on Deep Learning. *International Journal of Molecular Sciences*, 26(13), 6284. PMID: 40650062. https://doi.org/10.3390/ijms26136284

5. (Moqri, 2024) Moqri, M., Herzog, C., Poganik, J. R., Ying, K., Justice, J. N., et al. (2024). Validation of biomarkers of aging. *Nature Medicine*, 30, 360–372. PMID: 38355974. https://doi.org/10.1038/s41591-023-02784-9

6. (Levy, 2025) Levy, J. J., Diallo, A. B., Saldias Montivero, M. K., Gabbita, S., & Salas, L. A. (2025). Insights to aging prediction with AI based epigenetic clocks. *Epigenomics*, 17(1), 1–14. PMID: 39584810. https://doi.org/10.1080/17501911.2024.2432854

7. (Johnson, 2022) Johnson, A. A., English, B. W., Shokhirev, M. N., Sinclair, D. A., & Cuellar, T. L. (2022). Human age reversal: Fact or fiction? *Aging Cell*, 21(8), e13664. PMID: 35778957. https://doi.org/10.1111/acel.13664

8. (Richardson, 2025) Richardson, M., Brandt, C., Jain, N., Li, J. L., & Demanelis, K. (2025). Characterization of DNA methylation clock algorithms applied to diverse tissue types. *Aging*, 17(1). PMID: 39754638. https://doi.org/10.18632/aging.206182

9. (Herzog, 2025) Herzog, C. M. S., Redl, E., Barrett, J., Aminzadeh-Gohari, S., & Weber, D. D. (2025). Functionally enriched epigenetic clocks reveal tissue-specific discordant aging patterns in individuals with cancer. *Communications Medicine*, 5, 119. PMID: 40175686. https://doi.org/10.1038/s43856-025-00739-4

10. (Vetter, 2023) Vetter, V. M., Spieker, J., Sommerer, Y., Buchmann, N., & Kalies, C. H. (2023). DNA methylation age acceleration is associated with risk of diabetes complications. *Communications Medicine*, 3, 16. PMID: 36765171. https://doi.org/10.1038/s43856-023-00250-8

11. (Rutledge, 2022) Rutledge, J., Oh, H., & Wyss-Coray, T. (2022). Measuring biological age using omics data. *Nature Reviews Genetics*, 23, 715–727. https://doi.org/10.1038/s41576-022-00511-7

12. (Davydova, 2024) Davydova, E., Perenkov, A., & Vedunova, M. (2024). Building Minimized Epigenetic Clock by iPlex MassARRAY Platform. *Genes*, 15(4), 425. PMID: 38674360. https://doi.org/10.3390/genes15040425

13. (Shokhirev, 2025) Shokhirev, M. N., & Johnson, A. A. (2025). Using buccal methylomic data to create explainable aging clocks as well as classifiers and regressors for lifestyle and demographic factors. *Frontiers in Genetics*, 16, 1637186. PMID: 41104118. https://doi.org/10.3389/fgene.2025.1637186

14. (Tian, 2023) Tian, Y., Cropley, V., Maier, A. B., Lautenschlager, N. T., Breakspear, M., & Zalesky, A. (2023). Heterogeneous aging across multiple organ systems and prediction of chronic disease and mortality. *Nature Medicine*, 29, 1224–1232. https://doi.org/10.1038/s41591-023-02296-6

## Search Strategy

**Keywords used:**
1. "epigenetic clock DNA methylation biological age estimation deep learning" → PubMed
2. "Horvath clock GrimAge aging acceleration biomarker validation" → PubMed
3. "tissue-specific epigenetic clock DNA methylation aging organ" → PubMed
4. "exercise diet lifestyle intervention epigenetic clock age acceleration" → PubMed
5. "longevity cohort DNA methylation aging biomarker centenarian" → PubMed
6. "AltumAge deep neural network pan-tissue epigenetic clock" → PubMed
7. "neural network epigenetic clock tissue-specific DNA methylation aging" → OpenAlex (off-topic results)
8. "deep learning epigenetic aging clock" → ArXiv (network timeout)

**Date range:** 2013–2026 (majority from 2020–2026)
**Inclusion criteria:** English language; epigenetic/DNA methylation aging; human studies; methodological development
**Exclusion criteria:** Animal-only studies without translational relevance; non-aging outcomes

## Research Gaps Identified

1. Most existing clocks are trained on **blood tissue** only; cross-tissue generalisation remains poor
2. **Interventional studies** show modest effect sizes; detecting subtle lifestyle effects requires larger cohorts
3. **Deep learning models** achieve high correlation but suffer from poor extrapolation to extreme ages (longevity cohorts)
4. No consensus on **multi-tissue training strategy** for pan-tissue clocks
5. **Explainability** of neural network clocks is still limited despite recent XAI attempts (e.g., XAI-AGE)
