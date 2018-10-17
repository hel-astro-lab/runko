#!/bin/bash -l
#SBATCH -A 2017-12-57
#SBATCH -J pytests
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH -t 00-00:00:30
#SBATCH --nodes=1 
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1

module load git
module load anaconda
module load fftw
module load cdt/17.10
module load hdf5


cd $PLASMABOXDIR
aprun -n 1 -N 1 -d 1 python -m unittest discover -s tests/ -v

