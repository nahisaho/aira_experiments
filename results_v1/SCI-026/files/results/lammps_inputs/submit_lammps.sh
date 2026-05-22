#!/bin/bash
#SBATCH --job-name=LAMMPS_LPS_LCO
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00

module load lammps/2023.08.02

# Loop over temperatures for Arrhenius analysis
for T in 600 700 800 900 1000; do
  mkdir T${T}K && cd T${T}K
  sed "s/T_target equal 700/T_target equal ${T}/" ../lammps_md_interdiffusion.in > run.in
  mpirun -np $SLURM_NTASKS lmp < run.in > md_${T}K.log 2>&1
  cd ..
done
