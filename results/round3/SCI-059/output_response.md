Good, I have enough literature data. Now let me implement the full geostatistical framework.Now let me write the main experiment runner:Experiment complete. Let me check the figures and then write the deliverables:Both pass (report ≥ 850, paper ≥ 1500). Now write minimal tests and clean up:All 6 tests pass. Now clean up:All complete. Here is the summary:

---

## Experiment Complete: Geostatistical Disease Risk Framework

### Key Scientific Findings

1. **Strong spatial clustering of malaria prevalence** — Global Moran's I = **0.490** (z = 10.36, p < 0.001, permutation p < 0.002), confirming the necessity of spatial models over standard regression.

2. **Matérn variogram: range ≈ 0.776** — Spatial correlation persists over ~27% of the study domain. Nugget/Sill ratio = 0.24, indicating 76% of variance is spatially structured.

3. **Bayesian GP/SPDE outperforms Ordinary Kriging** — RMSE: 0.1261 vs 0.1296 (+2.7%); R²: 0.378 vs 0.343 (+10%), driven by automatic hyperparameter optimisation via marginal likelihood.

4. **SAR spatial lag coefficient ρ = 0.519** — ~52% of malaria risk explained by neighbouring sites, confirming strong spillover effects. SAR residual Moran's I = 0.197 (p < 0.001) indicates remaining nonlinear spatial structure.

5. **Spatiotemporal spline 5-fold CV**: RMSE = **0.328 ± 0.007**, R² = 0.064 ± 0.012 on dengue log-counts (Poisson stochasticity limits Gaussian model R²).

### Most Important Figure
![Malaria Risk Map — GP/SPDE predictions with uncertainty](figures/malaria_risk_map.png)

### Deliverables

**Source code:** 6 modules, 1,554 lines total  
**Tests:** 6/6 passing  
**Figures:** 8 (LGCP, Moran scatter, permutation test, variogram, 2× risk maps, spatiotemporal, comparison)  
**report.md:** 1,307 words | **paper.md:** 3,483 words (16 references, ≥30% post-2020)

### Limitations
- GP/SPDE uses Gaussian likelihood; Poisson/NegBin models needed for true count data
- Stationary covariance assumed; real landscapes require non-stationary extensions
- Results validated on synthetic data only; Malaria Atlas Project data validation is future work