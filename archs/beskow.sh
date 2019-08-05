#modules
module load cdt/19.06
#module load intel/18.0.0.128
#module switch PrgEnv-cray/5.2.82 PrgEnv-intel
module load cray-fftw/3.3.8.3
module load cray-hdf5/1.10.5.0
#module load anaconda/py36/4.3
#module switch gcc/4.8.1 gcc/7.3.0
#module load gcc/7.3.0
#module load cmake/3.15.5
module load craype-hugepages2M
module load mpi4py/3.0.2/py37
#module load git

#Switching to atp/2.1.3.
#Switching to cce/9.0.0.
#Switching to cray-fftw/3.3.8.3.
#Switching to cray-hdf5/1.10.5.0.
#Switching to cray-libsci/19.06.1.
#Switching to cray-mpich/7.7.8.
#Switching to craype/2.6.0.
#Switching to gcc/8.3.0.
#Switching to modules/3.2.11.2.
#Switching to perftools-base/7.1.0.
#Switching to pmi/5.0.14

# cray environment
export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

# enable dynamic library building with cray
export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

# python for CMAKE 
export PYTHON_EXECUTABLE=python

#prevent lockfiles
export HDF5_USE_FILE_LOCKING=FALSE

#no python .pyc files
export PYTHONDONTWRITEBYTECODE=true

# pybind installation
export LUSTREDIR=/cfs/klemming/nobackup/j/jnattila

export PYTHONDONTWRITEBYTECODE=true

unset PYTHONPATH
#eval "$(/cfs/klemming/nobackup/j/jnattila/miniconda3/bin/conda shell.bash hook)"
#export PATH=$PATH:/cfs/klemming/nobackup/j/jnattila/pkgs/lib/

#add libary paths
export RUNKODIR=$LUSTREDIR/runko
export PYTHONPATH=$PYTHONPATH:$RUNKODIR
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/python
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/analysis


