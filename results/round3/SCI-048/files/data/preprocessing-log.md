# Data Preprocessing Log

## Synthetic Data Generation

All data used in this study is synthetically generated.

### Burgers Equation Reference Solution
- **Method**: Method-of-lines with BDF integrator (scipy.integrate.solve_ivp)
- **Grid**: nx=256 points on x∈[-1,1], t_eval=101 points on [0,0.4]
- **Boundary conditions**: Dirichlet u(±1,t)=0, IC u(x,0)=-sin(πx)
- **Evaluation**: RegularGridInterpolator (bilinear)
- **Validity range**: t∈[0,0.4] (near-shock formation at t≈0.1/π²)

### Training Collocation Points
- **Method**: Latin Hypercube Sampling (LHS) via numpy.random.Generator
- **IC points**: 50 uniform x-values at t=0
- **BC points**: 20 Chebyshev-spaced t-values at x=±1
- **Collocation**: 400–1000 LHS points

### Operator Learning Dataset (Darcy Flow)
- **Input functions**: 250 samples = sum of 4 sinusoids with random frequencies (1–5) and amplitudes
- **Solutions**: Exact FEM solution via sparse tridiagonal solve
- **Split**: 200 train / 50 test
- **Grid**: nx=64 uniform on [0,1]

### Navier-Stokes Dataset
- **Reference**: Analytical Taylor-Green vortex (exact solution)
- **Training**: 400 random collocation + 100 IC points
- **Domain**: [0,2π]² × [0,0.5], Re=100

## No External Data Used
No external datasets or proprietary data were used in this study.
