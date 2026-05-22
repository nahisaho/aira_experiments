# Synthetic Data Preprocessing Log

## Global settings
- Random seed: `numpy.random.seed(42)` and `random.seed(42)` in every Python module.
- Figure export: PNG, 300 DPI, English-only labels, colorblind-friendly palette (`viridis` or `#0072B2`, `#E69F00`, `#009E73`, `#CC79A7`).
- Repository-relative output directories created automatically by each script: `results/` and `figures/`.

## 1. Gait Parkinson module
- Samples: 200 subjects (100 Parkinson disease, 100 healthy control).
- Synthetic signal design: 12-second accelerometer and gyroscope traces at 50 Hz.
- Parkinson group shifts: longer stride interval, higher stride variability, stronger 3-8 Hz freezing component, lower step regularity, slower turning peaks, lower cadence, lower walking speed.
- Extracted features:
  - stride length variability as coefficient of variation of stride intervals
  - gait asymmetry index from left/right stride intervals
  - freeze of gait power ratio using FFT energy in 3-8 Hz divided by 0.5-3 Hz
  - step regularity by autocorrelation peak
  - turning speed from gyroscope peak detection
  - cadence and walking speed estimates
- Modeling: RandomForestClassifier with 5-fold stratified cross-validation.

## 2. Voice ALS module
- Samples: 150 recordings (50 healthy, 50 early ALS, 50 moderate ALS).
- Synthetic acoustic feature priors vary by severity for jitter, shimmer, HNR, F0, F0 standard deviation, speaking rate, and ALSFRS-R proxy.
- MFCC generation: 13 mean coefficients and 13 standard deviations per sample with class-dependent shifts.
- Preprocessing: standardization inside SVM and SVR pipelines.
- Tasks: multi-class classification and continuous regression of ALSFRS-R proxy.

## 3. Touchscreen cognition module
- Samples: 180 sessions (60 healthy, 60 MCI, 60 mild dementia).
- Features generated from class-dependent Gaussian priors:
  - inter-tap interval mean and CV
  - tap accuracy distance
  - swipe velocity and acceleration
  - long-press duration
  - double-tap timing variability
  - error rate and correction frequency
  - typing rhythm entropy
- Modeling: XGBoost if installed, otherwise GradientBoostingClassifier; 5-fold stratified cross-validation.

## 4. Change-point detection module
- Patients: 3 synthetic longitudinal trajectories across 24 months.
- Ground truth:
  - Patient 1 stable (no change point)
  - Patient 2 change at month 12
  - Patient 3 change at month 8
- Composite biomarker score: baseline + mild trend + Gaussian noise + optional post-change offset.
- Detection methods:
  - CUSUM with early baseline calibration
  - PELT using `ruptures` if available, otherwise a custom two-segment SSE optimization
  - simplified BOCPD using Gaussian z-score evidence
  - ensemble confirmation when 2 of 3 methods agree within ±1 month

## 5. Multimodal fusion module
- Assumption logged: request describes a 4-class problem; implementation uses 40 samples for each of Healthy, Parkinson, ALS, and Cognitive groups (160 total).
- Modality blocks:
  - gait: 7 features inspired by gait analysis outputs
  - voice: 8 features inspired by dysarthria monitoring outputs
  - touchscreen: 10 features inspired by cognitive interaction outputs
- Fusion strategies:
  - early fusion via concatenation + Random Forest
  - late fusion via weighted averaging of modality-specific probabilities
  - attention proxy via MLP on concatenated features
- Composite NDD-Score: scaled 0-100 from disease probability plus severity proxy.

## 6. Validation strategy module
- Patients: 120 matched synthetic cases (40 Parkinson, 40 ALS, 40 cognitive impairment).
- Clinical endpoint simulation:
  - digital gait score vs UPDRS-III
  - digital voice score vs ALSFRS-R impairment transform
  - digital touchscreen score vs MoCA impairment transform
- Bootstrap: 1,000 resamples for Pearson correlation confidence intervals.
- Agreement: concordance correlation coefficient and Bland-Altman analysis on aligned score scales.
- Sensitivity to change: SRM from simulated baseline and follow-up scores.
- Power analysis: normal approximation for two-arm trial detecting 20% score change at 80% power and alpha 0.05.
