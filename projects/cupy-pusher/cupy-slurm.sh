#!/bin/bash
#SBATCH --account=project_462001058
##SBATCH --partition=standard-g
#SBATCH --partition=dev-g
#
#SBATCH --job-name=cupy
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --open-mode=truncate
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-00:10:00       # Run time (d-hh:mm:ss)

# Loads correct modules and sets up PYTHONPATH.
source runko/archs/lumi-load-runko-env
module load CuPy/13.4.1-cpeGNU-24.03-rocm-6.2.2.

export OMP_NUM_THREADS=8

# Make output directory
OUTDIR="runs/${SLURM_JOB_NAME}-$SLURM_JOB_ID"
mkdir $OUTDIR
# Go
srun python ./runko/projects/cupy-pusher/pic-cupy.py $OUTDIR

# Save output files in handier place
mv ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out $OUTDIR/slurm.out
mv ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err $OUTDIR/slurm.err