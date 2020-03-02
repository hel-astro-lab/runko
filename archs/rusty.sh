#modules
module load slurm/18.08.8
module load gcc/7.4.0
module load openmpi2/2.1.6-hfi
module load python3/3.7.3
module load python3-mpi4py/3.7.3-openmpi2
module load cmake/3.14.5
module load modules-nix/20200227-37
module load nix/ffmpeg/3.4.7
module load nix/fftw/3.3.8
module load lib/hdf5/1.10.5-openmpi2


# cray environment
export CXX=mpic++
export CC=mpicc

# python for CMAKE 
export PYTHON_EXECUTABLE=python

#prevent lockfiles
export HDF5_USE_FILE_LOCKING=FALSE

#no python .pyc files
export PYTHONDONTWRITEBYTECODE=true

# pybind installation
export LUSTREDIR=/mnt/home/$USER/ceph

#unset PYTHONPATH
#eval "$(/cfs/klemming/nobackup/j/jnattila/miniconda3/bin/conda shell.bash hook)"
#export PATH=$PATH:/cfs/klemming/nobackup/j/jnattila/pkgs/lib/

#add libary paths
export RUNKODIR=$LUSTREDIR/runko
export PYTHONPATH=$PYTHONPATH:$RUNKODIR
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/bindings/old

