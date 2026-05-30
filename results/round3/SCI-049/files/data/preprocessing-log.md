# Preprocessing Log

Synthetic data were generated with deterministic seeds, standardized inside each detector, and streamed in windows of 100 samples. Training used the first 60% of samples; the last 40% were held out for evaluation. When Page-Hinkley drift detection fired, models were retrained on the lowest-score 80% of the most recent 800 samples to reduce contamination from anomalies.
