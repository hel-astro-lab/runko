#!/bin/bash
#SBATCH -A 2017-12-57
#SBATCH -J pytests
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 00-00:00:30
#SBATCH --nodes=1 
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1

cd $PLASMABOXDIR
aprun -n 1 -N 1 -d 1 python3 -m unittest discover -s tests/ -v

