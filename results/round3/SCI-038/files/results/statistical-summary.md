# Statistical Summary

## Assumption Checks
For the paired seed study comparing nearest-neighbor and GA+2opt mission delta-v, the paired differences passed the Shapiro-Wilk normality check (p = 0.844), so a paired t-test was used. For the capture-method comparison between harpoon and robotic arm across the six tested spin-rate conditions, paired differences also passed the Shapiro-Wilk check (p = 0.476).

## Test 1: Sequence Optimization Across Seeds
- Nearest-neighbor: 14.834 ± 4.614 km/s; 95% CI [9.992, 19.677]
- GA+2opt: 14.807 ± 4.966 km/s; 95% CI [9.596, 20.018]
- Mean paired reduction (NN - GA): 0.028 ± 0.444 km/s; 95% CI [-0.438, 0.493]
- Paired t-test: t = 0.152, p = 0.8850
- Effect size: Cohen's dz = 0.062

Interpretation: the baseline run showed a small improvement for GA+2opt, but the six-seed study indicates that the advantage is not statistically reliable for this simplified synthetic setting. In practical terms, the greedy baseline is already strong, and the evolutionary refinement mainly provides a mechanism for occasional local improvements rather than a guaranteed mission-level gain.

## Test 2: Capture Method Difference
- Mean harpoon minus robotic-arm success probability: 0.169 ± 0.157; 95% CI [0.004, 0.334]
- Paired t-test: t = 2.636, raw p = 0.0462
- Bonferroni-corrected p (2 planned comparisons): 0.0924
- Effect size: Cohen's dz = 1.076

Interpretation: harpoon-like capture retains a materially higher modeled success probability than robotic-arm capture as rotation rate increases, but the corrected p-value shows that the evidence is suggestive rather than definitive under the present small-sample synthetic experiment. The effect size remains large, so the result is still operationally relevant for mechanism screening.
