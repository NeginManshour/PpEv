

echo "### Starting at: $(date) ###"

module load miniconda3
source activate haddock3

export TMPDIR=$SLURM_SUBMIT_DIR/${SLURM_JOB_ID}-${SLURM_JOB_NAME}
mkdir -p $TMPDIR
cd $TMPDIR

python /home/nmn5x/data/IOU/haddock3/Data_haddock/TF/7zx4_TF/run_haddock_score.py /home/nmn5x/data/IOU/haddock3/Data_haddock/TF/7zx4_TF/pdb/ 1000 5

echo "### Ending at: $(date) ###"
