# Preprocessing Log

Synthetic TMS-EEG responses were generated with fixed random seeding (`np.random.seed(42)`), channel-wise baseline normalization, Gaussian temporal smoothing, and additive colored plus heavy-tailed noise to emulate biologically plausible variability. DOC feature matrices were derived from out-of-fold simulations only; no fold-specific label statistics were leaked into feature construction or classifier fitting. Transfer-entropy heatmaps were computed on binarized Markov-chain state trajectories sampled from anesthesia depth–dependent TPMs.
