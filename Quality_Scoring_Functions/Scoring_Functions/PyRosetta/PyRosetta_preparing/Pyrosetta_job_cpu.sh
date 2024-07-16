#!/bin/bash
#SBATCH -p general
#SBATCH --mem 150G
#SBATCH -c 4
#SBATCH --time 02-00:00:00

#SBATCH --job-name=P_F_8FZM

echo "### Starting at: $(date) ###"

module load miniconda3
source activate pyrosetta

python D020_pose_scoring_bulk.py ./pdb

echo "### Ending at: $(date) ###"

