# Preprocessing Log

- Random seeds fixed with numpy default_rng(42).
- 100x100 grid and 50 acquisitions with 24-day repeat cycle.
- Synthetic signals: interseismic trend, slow slip event, annual/semi-annual deformation, atmospheric delays, and sensor noise.
- LOS geometries represent ascending and descending Sentinel-1-like tracks.
- GPS control points were sampled from the true ENU field with small Gaussian perturbations.
