#!/bin/bash
#SBATCH --account=project_462001058
##SBATCH --partition=standard-g
#SBATCH --partition=dev-g
#
#SBATCH --job-name=cherenkov
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --open-mode=truncate
#
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-00:15:00       # Run time (d-hh:mm:ss)

# Loads correct modules and sets up PYTHONPATH.
source runko/archs/lumi-load-runko-env

# Required to choose correct GPU for each task.
cat << EOF > select_gpu
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID

exec \$*
EOF

chmod +x ./select_gpu

CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

export OMP_NUM_THREADS=6
export MPICH_GPU_SUPPORT_ENABLED=1

# Make output directory
OUTDIR="runs/${SLURM_JOB_NAME}-$SLURM_JOB_ID"
mkdir $OUTDIR
# Go
srun --cpu-bind=${CPU_BIND} ./select_gpu python ./runko/projects/pic2-tube/pic-cherenkov.py $OUTDIR

rm -f ./select_gpu

# Save output files in handier place
mv ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out $OUTDIR/slurm.out
mv ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err $OUTDIR/slurm.err