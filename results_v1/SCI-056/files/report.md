# DRAFT — NOT FOR DISTRIBUTION

## COVID-19 Wave 6/7 retrospective case study

Timestamp: 2026-05-22T17:28:23.888941+00:00

## Methods

This analysis generated synthetic daily COVID-19 case curves for Japan's 6th wave (Omicron BA.1, January-March 2022) and 7th wave (BA.5, July-September 2022) using SEIR-calibrated epidemic shapes, realistic weekday reporting noise, age-specific case allocation (20%/30%/30%/20%), booster rollout trajectories, and lagged severity assumptions. Model fitting used least-squares estimation for SIR, SEIR, and age-structured SEIR models. Model comparison used AIC and BIC.

## Results

- Wave 6 fitted SEIR R0: 5.00
- Wave 7 fitted SEIR R0: 8.00
- Wave 6 combined-intervention peak reduction vs no intervention: 92.9%
- Wave 7 combined-intervention peak reduction vs no intervention: 84.9%
- Best model for Wave 6: Age-structured SEIR
- Best model for Wave 7: Age-structured SEIR

## Discussion

The synthetic case study reproduces the higher transmissibility of BA.5 relative to BA.1, the contribution of booster coverage to severity reduction, and the added value of combining contact reduction with vaccination. As requested, outputs are figure-ready datasets rather than rendered plots. The age-structured model improves severity alignment because severe outcomes are concentrated in older adults even when total case shares are more balanced.

## File inventory

- `results/covid_case_study.py`
- `results/covid_wave6_results.json`
- `results/covid_wave7_results.json`
- `results/scenario_comparison.json`
- `results/model_comparison.json`
- `results/statistical-summary.md`
- `data/preprocessing-log.md`
- `logs/process-log.jsonl`

## Key findings

```json
{
  "R0_estimates": {
    "wave6": 5.0,
    "wave7": 8.0
  },
  "intervention_effectiveness": {
    "wave6_peak_reduction_pct": 92.865156,
    "wave6_total_reduction_pct": 89.516907,
    "wave7_peak_reduction_pct": 84.862776,
    "wave7_total_reduction_pct": 69.4294,
    "wave6_no_intervention_total_cases": 47695845.342555,
    "wave7_no_intervention_total_cases": 39253400.932766
  },
  "model_comparison_results": {
    "wave6_best_model": "Age-structured SEIR",
    "wave7_best_model": "Age-structured SEIR"
  },
  "lessons_learned": [
    "Higher BA.5 transmissibility increased the fitted basic reproduction number from the BA.1 wave to the BA.5 wave.",
    "Vaccination alone reduced severe outcomes more strongly than it reduced transmission, consistent with booster-era Omicron experience.",
    "Combined interventions were required for Reff to fall sustainably below 1 in both waves.",
    "Age-structured modelling improved severity fit because hospitalization and death burden concentrated in older adults despite similar case shares."
  ]
}
```
