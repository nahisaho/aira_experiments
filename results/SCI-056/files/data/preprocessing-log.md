# Preprocessing log

- Analysis type: synthetic retrospective COVID-19 case study for Japan Wave 6 and Wave 7.
- Random seed: numpy default_rng(20240219).
- Data source: synthetic data calibrated to published-scale epidemic summaries (daily peaks, cumulative cases, booster rollout, age severity pattern).
- Transformation steps:
  1. Construct booster coverage and contact reduction curves with logistic functions.
  2. Calibrate SEIR-based epidemic curves to target peak size, total cases, and peak timing for each wave.
  3. Add realistic weekday reporting variation and mild log-normal observation noise.
  4. Allocate cases to four age bands using time-varying shares around 20/30/30/20.
  5. Apply lagged severity rates to obtain hospitalizations, ICU admissions, and deaths.
- Reproducibility: numpy/scipy only; no external datasets required.
