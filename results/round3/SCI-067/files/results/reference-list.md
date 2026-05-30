# Reference List

## Included Studies (Screened and Included)

1. (Köck, 2023) Köck, B., Friedl, A., Serna Loaiza, S., Wukovits, W., & Mihalyi-Schneider, B. (2023). Automation of Life Cycle Assessment—A Critical Review of Developments in the Field of Life Cycle Inventory Analysis. *Sustainability*, 15(6), 5531. https://doi.org/10.3390/su15065531

2. (Wernet, 2016) Wernet, G., Bauer, C., Steubing, B., Reinhard, J., Moreno-Ruiz, E., & Weidema, B. (2016). The ecoinvent database version 3 (part I): Overview and methodology. *The International Journal of Life Cycle Assessment*, 21(9), 1218–1230. https://doi.org/10.1007/s11367-016-1087-8

3. (Mutel, 2017) Mutel, C. (2017). Brightway2: an advanced framework for life cycle assessment in Python. *Journal of Open Source Software*, 2(12), 472. https://doi.org/10.21105/joss.00472

4. (Gutsch, 2023) Gutsch, M., & Leker, J. (2023). Costs, carbon footprint, and environmental impacts of lithium-ion batteries — From cathode active material synthesis to cell manufacturing and recycling. *Applied Energy*, 352, 122132. https://doi.org/10.1016/j.apenergy.2023.122132

5. (Llamas-Orozco, 2023) Llamas-Orozco, J. A., et al. (2023). Estimating the environmental impacts of global lithium-ion battery supply chain: A temporal, geographical, and technological perspective. *PNAS Nexus*, 2(11), pgad361. https://doi.org/10.1093/pnasnexus/pgad361

6. (Nguyen, 2023) Nguyen, Q., Diaz-Rainey, I., et al. (2023). Scope 3 emissions: Data quality and machine learning prediction accuracy. *PLOS Climate*, 2(11), e0000208. https://doi.org/10.1371/journal.pclm.0000208

7. (Jain, 2023) Jain, A., Padmanaban, M., et al. (2023). Supply chain emission estimation using large language models. *arXiv preprint*. https://doi.org/10.48550/arXiv.2308.01741

8. (Jain, 2024) Jain, A., Padmanaban, M., et al. (2024). A Framework for Emission Reduction in Scope 3 Climate Actions using Domain-adapted Foundation Model. *CODS-COMAD 2024*. https://doi.org/10.1145/3632410.3632465

9. (Serafeim, 2022) Serafeim, G., & Vélez Caicedo, G. (2022). Machine Learning Models for Prediction of Scope 3 Carbon Emissions. *Harvard Business School Working Paper*, No. 22-080. https://www.hbs.edu/ris/Publication%20Files/22%20080_035d70d9-3acf-4faa-aa93-534e52a52d0e.pdf

10. (Lai, 2022) Lai, X., et al. (2022). Life Cycle Assessment of Lithium-ion Batteries for Carbon-peaking and Carbon-neutrality. *Journal of Mechanical Engineering*, 58(22), 3–20. https://doi.org/10.3901/JME.2022.22.003

11. (Huijbregts, 2017) Huijbregts, M. A. J., et al. (2017). ReCiPe2016: a harmonised life cycle impact assessment method at midpoint and endpoint level. *The International Journal of Life Cycle Assessment*, 22, 138–147. https://doi.org/10.1007/s11367-016-1246-y

12. (Saltelli, 2020) Saltelli, A., et al. (2020). Five ways to ensure that models serve society: a manifesto. *Nature*, 582, 482–484. https://doi.org/10.1038/s41586-020-2484-8

## Research Gaps Identified

1. **Automated process tree extraction**: Most LCA automation remains at the data-linking stage; end-to-end NLP-driven tree construction from unstructured product documents is understudied.
2. **Integrated uncertainty + hotspot pipeline**: Monte Carlo and Taylor expansion methods exist separately; automated joint pipelines with real-time hotspot flagging are rare.
3. **Scope 3 emission estimation via ML in LCA context**: Existing ML methods focus on financial proxies, not LCA-process-tree-based estimation.
4. **EV battery with full supply-chain LCA automation**: Most automated approaches cover single-stage analysis; multi-tier Scope 3 integration is limited.
5. **Ecoinvent automated semantic matching**: TF-IDF/cosine similarity matching exists but precision on ambiguous process names is low (~70–80%).
