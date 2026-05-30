# Experimental Design

Objective: evaluate automated quality control for CERN/LIGO-scale scientific streams using synthetic multivariate detector-like data with controllable anomalies, structural breaks, and concept drift.

## Dataset design
- 5,000 time steps and 8 channels representing detector subsystems.
- Point anomalies: 1.0% of timestamps.
- Contextual anomalies: 0.5% of timestamps.
- Changepoints: three locations near 24%, 48%, and 70% of the stream.
- Concept drift: gradual regime shift starting at t=3500.
- Physical constraints: positivity on key channels, additive consistency, derived-channel consistency, coupled-channel consistency, and aggregate-balance consistency.
- Noise model: Gaussian noise calibrated to approximately 20 dB SNR; realism guard reruns with more noise if any metric becomes effectively perfect.

## Methods
1. PELT changepoint detection using `ruptures`.
2. Isolation Forest for pointwise unsupervised anomaly scoring.
3. Simplified Deep SVDD with a lightweight PyTorch autoencoder and center loss.
4. Physics-based constraint scoring from residual violations.
5. Page-Hinkley drift detection on streaming batch statistics.
6. Weighted ensemble fusion of normalized detector outputs.
7. SHAP-based explanation of flagged anomalies using Kernel SHAP with deterministic fallback.

## Streaming protocol
- Initial training window: first 60% of samples.
- Evaluation window: last 40% of samples.
- Processing batch size: 100 samples.
- Drift adaptation: retrain on the lowest-score 80% of the most recent 800 samples after a Page-Hinkley alarm.

## Metrics
- Overall detection: F1, precision, recall, AUROC.
- Changepoints: mean detection delay and false alarm rate.
- Drift: detection delay and number of retraining triggers.
- Stability: 5-fold seeded cross-validation with mean ± standard deviation.
