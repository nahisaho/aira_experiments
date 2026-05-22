# Simulation Protocol Summary

## Phase 1: System Preparation
| Step | Tool | Duration | Notes |
|------|------|----------|-------|
| Topology generation | GROMACS `pdb2gmx` / moltemplate | — | OPLS-AA + scaled charges |
| Box setup | `gmx insert-molecules` / `packmol` | — | Random placement |
| Energy minimization | `gmx mdrun` (steep) | 50,000 steps | Fmax < 100 kJ/mol/nm |
| NVT equilibration | `gmx mdrun` | 500 ps | T = 298.15 K, v-rescale |
| NPT equilibration | `gmx mdrun` | 2 ns | T = 298.15 K, P = 1 bar |

## Phase 2: Production Runs
| Property | Ensemble | Duration | Δt | Save freq |
|----------|----------|----------|----|-----------|
| Structural (RDF, CN) | NPT | 20 ns | 2 fs | 1 ps |
| KB integrals | NPT | 50 ns | 2 fs | 1 ps |
| Diffusion (MSD) | NVE | 20 ns | 1 fs | 0.1 ps |
| Conductivity (GK) | NVE | 50 ns | 1 fs | 0.1 ps |
| Solvation FE (TI) | NPT | 5 ns × 21 λ | 2 fs | 1 ps |

## Phase 3: Analysis
| Analysis | Method | Script |
|----------|--------|--------|
| Activity coefficient | Kirkwood-Buff integrals | `scripts/01_kirkwood_buff.py` |
| Diffusion coefficient | Einstein relation (MSD) | `scripts/02_transport_msd.py` |
| Conductivity | Green-Kubo (current ACF) | `scripts/03_green_kubo_conductivity.py` |
| Solvation structure | RDF + coordination number | `scripts/04_solvation_structure.py` |
| Solvation free energy | Thermodynamic integration | `scripts/05_solvation_free_energy.py` |
| Anomalous transport | Fractional Stokes-Einstein | `scripts/06_anomalous_transport.py` |

## Computational Resources (Estimated)
- Single system (20 ns NPT): ~24 CPU-hours on 64 cores
- Green-Kubo (50 ns NVE): ~80 CPU-hours on 64 cores  
- TI (21 windows × 5 ns): ~120 CPU-hours on 64 cores
- Total for 5 systems: ~1,500–2,000 CPU-hours
