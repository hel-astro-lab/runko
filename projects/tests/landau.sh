#!/bin/bash
#SBATCH -A SNIC2018-5-16
#SBATCH -J lan2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 01-02:00:00
#SBATCH -n 1
#SBATCH -c 6

# export module libraries
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/plasma-toolbox
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/plasma-toolbox/python

# activate threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# go to working directory
cd /pfs/nobackup/home/n/natj/plasma-toolbox/projects/tests/

mpirun -np 1 python plasma.py --conf config-landau.ini
