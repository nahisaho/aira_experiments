# Preprocessing Log

- Synthetic dataset generator: `sklearn.datasets.make_classification`.
- Quantum-kernel benchmark: 200 samples, 4 features, additive Gaussian noise with standard deviation 0.1, stratified 5-fold CV.
- Encoding benchmark: 160 samples, 4 features, additive Gaussian noise with standard deviation 0.1, stratified 5-fold CV.
- Feature scaling:
  - Quantum feature maps: MinMax scaling to [-1, 1].
  - Classical RBF-SVM: Standard scaling.
- Random seeds fixed in all experiment modules for reproducibility.
- No real-world data, patient data, or proprietary data used.
