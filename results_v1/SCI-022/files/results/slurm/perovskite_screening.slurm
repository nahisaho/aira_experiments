#!/bin/bash
#SBATCH --job-name=perovskite_screening
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

module load intel/2023a VASP/6.4.1
conda activate perovskite

# AiiDA daemon must be running
verdi daemon start 4

# Submit workflow
python -c "
from perovskite_screener.workflow import build_screening_workflow
from perovskite_screener.materials_database import get_all_candidates
import json

candidates = get_all_candidates()
wf = build_screening_workflow(candidates)
wf.save_json('results/workflow_definition.json')
print(f'Workflow {wf.wf_id} with {len(wf.fireworks)} fireworks saved.')
"
