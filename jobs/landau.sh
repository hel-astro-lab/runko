#!/bin/bash
#SBATCH -A SNIC2018-5-16
#SBATCH -J landau
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 00-00:15:00
#SBATCH -n 1
#SBATCH -c 4

# export module libraries
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/plasma-toolbox
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/plasma-toolbox/python

# activate threading
export OMP_NUM_THEADS=$SLURM_CPUS_PER_TASK

# go to working directory
cd /pfs/nobackup/home/n/natj/plasma-toolbox/projects/tests/

mpirun -np 1 python plasma.py
#srun python plasma.py
