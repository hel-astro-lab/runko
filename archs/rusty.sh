#modules
module load gcc/11.4.0
module load cmake
module load openmpi/4.0.7
module load python-mpi/3.10.10
module load hdf5/mpi-1.10.9
module load texlive


# cray environment
export CXX=mpic++
export CC=mpicc

# python for CMAKE 
export PYTHON_EXECUTABLE=python3

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
export RUNKODIR=$LUSTREDIR/runko-gpu
export PYTHONPATH=$PYTHONPATH:$RUNKODIR
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/bindings/old

