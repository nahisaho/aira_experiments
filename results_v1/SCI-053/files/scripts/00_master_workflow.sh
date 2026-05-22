#!/bin/bash
# =============================================================
# Master Workflow: Concentrated Electrolyte MD Simulation
# =============================================================
# System: LiPF6 in EC:DMC (1:1 vol)
# Force Field: OPLS-AA with ECC (scaled charges, q_scale=0.8)
# Software: GROMACS 2023+
#
# Usage: bash scripts/00_master_workflow.sh [system_name]
# =============================================================

set -euo pipefail

SYSTEM=${1:-"ec_dmc_lipf6_1M"}
NPROC=${2:-8}
GMX="gmx"

echo "============================================================"
echo "Concentrated Electrolyte MD Simulation Pipeline"
echo "System: ${SYSTEM}"
echo "Processors: ${NPROC}"
echo "============================================================"

# ---------------------------------------------------------------
# Phase 0: System Setup
# ---------------------------------------------------------------
echo ""
echo "=== Phase 0: System Setup ==="

# 0.1 Generate molecular topologies
# (Assumes OPLS-AA itp files exist for EC, DMC, LiPF6)
# ${GMX} pdb2gmx -f ${SYSTEM}.pdb -o ${SYSTEM}_processed.gro \
#     -p topol.top -ff oplsaa

# 0.2 Insert molecules with packmol
# packmol < packmol_input.inp

# 0.3 Create box and add molecules
# ${GMX} editconf -f ${SYSTEM}.gro -o ${SYSTEM}_box.gro \
#     -c -d 1.0 -bt cubic

echo "  [NOTE] System setup requires topology files."
echo "  Place .itp and .gro files in the working directory."

# ---------------------------------------------------------------
# Phase 1: Energy Minimization
# ---------------------------------------------------------------
echo ""
echo "=== Phase 1: Energy Minimization ==="

# ${GMX} grompp -f results/gromacs_mdp/em.mdp \
#     -c ${SYSTEM}_box.gro -p topol.top -o em.tpr -maxwarn 2
# ${GMX} mdrun -v -deffnm em -ntomp ${NPROC}

echo "  [TEMPLATE] Use em.mdp with steep integrator, 50000 steps"

# ---------------------------------------------------------------
# Phase 2: NVT Equilibration (500 ps)
# ---------------------------------------------------------------
echo ""
echo "=== Phase 2: NVT Equilibration ==="

# ${GMX} grompp -f results/gromacs_mdp/nvt_equil.mdp \
#     -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
# ${GMX} mdrun -deffnm nvt -ntomp ${NPROC}

echo "  [TEMPLATE] V-rescale thermostat, T=298.15K, 500 ps"

# ---------------------------------------------------------------
# Phase 3: NPT Equilibration (2 ns)
# ---------------------------------------------------------------
echo ""
echo "=== Phase 3: NPT Equilibration ==="

# ${GMX} grompp -f results/gromacs_mdp/npt_production.mdp \
#     -c nvt.gro -r nvt.gro -p topol.top -o npt_equil.tpr \
#     -maxwarn 2
# # Override nsteps for 2 ns equilibration
# ${GMX} mdrun -deffnm npt_equil -ntomp ${NPROC} -nsteps 1000000

echo "  [TEMPLATE] Parrinello-Rahman, P=1 bar, 2 ns"

# ---------------------------------------------------------------
# Phase 4: NPT Production (20 ns) - Structural Properties
# ---------------------------------------------------------------
echo ""
echo "=== Phase 4: NPT Production (Structure) ==="

# ${GMX} grompp -f results/gromacs_mdp/npt_production.mdp \
#     -c npt_equil.gro -p topol.top -o npt_prod.tpr -maxwarn 2
# ${GMX} mdrun -deffnm npt_prod -ntomp ${NPROC}

echo "  [TEMPLATE] 20 ns NPT, save every 1 ps"

# Structural analysis
echo "  Running structural analysis..."
# ${GMX} rdf -f npt_prod.xtc -s npt_prod.tpr -o rdf_all.xvg
# ${GMX} density -f npt_prod.xtc -s npt_prod.tpr -o density.xvg

# ---------------------------------------------------------------
# Phase 5: NVE Production (50 ns) - Transport Properties
# ---------------------------------------------------------------
echo ""
echo "=== Phase 5: NVE Production (Transport) ==="

# ${GMX} grompp -f results/gromacs_mdp/nve_transport.mdp \
#     -c npt_prod.gro -p topol.top -o nve_prod.tpr -maxwarn 2
# ${GMX} mdrun -deffnm nve_prod -ntomp ${NPROC}

echo "  [TEMPLATE] 50 ns NVE, 1 fs timestep, velocities every 0.1 ps"

# ---------------------------------------------------------------
# Phase 6: Analysis Pipeline
# ---------------------------------------------------------------
echo ""
echo "=== Phase 6: Analysis ==="

echo "  Running Kirkwood-Buff analysis..."
python3 scripts/01_kirkwood_buff.py

echo "  Running MSD diffusion analysis..."
python3 scripts/02_transport_msd.py

echo "  Running Green-Kubo conductivity..."
python3 scripts/03_green_kubo_conductivity.py

echo "  Running solvation structure analysis..."
python3 scripts/04_solvation_structure.py

echo "  Running solvation free energy analysis..."
python3 scripts/05_solvation_free_energy.py

echo "  Running anomalous transport analysis..."
python3 scripts/06_anomalous_transport.py

# ---------------------------------------------------------------
# Phase 7: Summary
# ---------------------------------------------------------------
echo ""
echo "============================================================"
echo "Pipeline Complete!"
echo "============================================================"
echo ""
echo "Output files:"
echo "  results/kb_analysis_results.json"
echo "  results/diffusion_results.json"
echo "  results/green_kubo_conductivity.json"
echo "  results/solvation_analysis.json"
echo "  results/solvation_free_energy.json"
echo "  results/anomalous_transport.json"
echo ""
echo "Figures:"
echo "  figures/kirkwood_buff_analysis.png"
echo "  figures/diffusion_analysis.png"
echo "  figures/green_kubo_conductivity.png"
echo "  figures/solvation_structure.png"
echo "  figures/anomalous_transport.png"
echo ""
echo "See report.md for full analysis summary."
