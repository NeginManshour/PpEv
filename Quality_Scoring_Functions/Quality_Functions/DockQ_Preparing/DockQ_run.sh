#!/bin/bash
#SBATCH -p gpu4
#SBATCH --mem 50G
#SBATCH -c 8
#SBATCH --gres gpu:2
#SBATCH --time 00-01:00:00

echo "### Starting at: $(date) ###"

module load miniconda3
source activate Evaluation

python run_dockq_bulk.py ./TB/template_based ./8CCW_AB.pdb   # As a sample : template_based model of 8CCW with chains A and B

echo "### Ending at: $(date) ###"

