#!/bin/bash
# Thermodynamic Integration Workflow for Li+ Solvation Free Energy
# System: Li+ in EC:DMC (1:1 vol)
# Protocol: Separate electrostatic and vdW decoupling

set -e

SYSTEM="ec_dmc_lipf6"
N_LAMBDA=21
T=298.15

echo "============================================"
echo "TI Solvation Free Energy Calculation"
echo "============================================"

# Step 1: Prepare topology with Li+ as separate molecule type
echo "Step 1: Preparing system..."
# gmx pdb2gmx -f ${SYSTEM}.pdb -o ${SYSTEM}.gro -p topol.top

# Step 2: Equilibrate reference state (lambda=0)
echo "Step 2: Equilibrating reference state..."
# gmx grompp -f em.mdp -c ${SYSTEM}.gro -p topol.top -o em.tpr
# gmx mdrun -deffnm em
# gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
# gmx mdrun -deffnm nvt
# gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr
# gmx mdrun -deffnm npt

# Step 3: Run electrostatic decoupling windows
echo "Step 3: Running electrostatic decoupling..."
for i in $(seq 0 $((N_LAMBDA-1))); do
    lambda=$(echo "scale=4; $i / ($N_LAMBDA - 1)" | bc)
    dir="elec/lambda_${i}"
    mkdir -p $dir
    echo "  Window $i: lambda = $lambda"
    # gmx grompp -f ${dir}/ti.mdp -c npt.gro -p topol.top -o ${dir}/ti.tpr
    # gmx mdrun -deffnm ${dir}/ti -ntomp 4
done

# Step 4: Run vdW decoupling windows
echo "Step 4: Running vdW decoupling..."
for i in $(seq 0 $((N_LAMBDA-1))); do
    lambda=$(echo "scale=4; $i / ($N_LAMBDA - 1)" | bc)
    dir="vdw/lambda_${i}"
    mkdir -p $dir
    echo "  Window $i: lambda = $lambda"
    # gmx grompp -f ${dir}/ti.mdp -c npt.gro -p topol.top -o ${dir}/ti.tpr
    # gmx mdrun -deffnm ${dir}/ti -ntomp 4
done

# Step 5: Analyze with gmx bar
echo "Step 5: Analyzing free energies..."
# gmx bar -f elec/lambda_*/ti.xvg -o elec_bar.xvg
# gmx bar -f vdw/lambda_*/ti.xvg -o vdw_bar.xvg

echo "Done!"
echo "Use scripts/05_solvation_free_energy.py for final analysis"
