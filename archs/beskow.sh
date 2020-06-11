#modules
module load cdt/19.06
module load intel/18.0.0.128
module swap PrgEnv-cray PrgEnv-intel
module load cray-fftw/3.3.8.3
module load mpi4py/3.0.2/py37


# cray environment
export CXX=CC
export CC="cc -D_Float128=__float128"
export FC=ftn
export MPICC="cc -D_Float128=__float128"
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

#unset PYTHONPATH
#eval "$(/cfs/klemming/nobackup/j/jnattila/miniconda3/bin/conda shell.bash hook)"
#export PATH=$PATH:/cfs/klemming/nobackup/j/jnattila/pkgs/lib/

#add libary paths
export RUNKODIR=$LUSTREDIR/runko
export PYTHONPATH=$PYTHONPATH:$RUNKODIR
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/bindings/old
