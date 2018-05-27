#!/bin/bash
#SBATCH -A SNIC2018-5-16
#SBATCH -J pytests
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 00-00:00:30
#SBATCH -n 1
#SBATCH -c 2
srun python -m unittest discover -s tests/ -v
