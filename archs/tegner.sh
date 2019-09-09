module load cmake/3.9.6
module load git
module load gcc/6.2.0
module load hdf5/1.10.0-gcc-6.2
module load anaconda/py35/4.2.0

#module unload fftw/3.3.7-gcc-7.2-openmpi-3.0-single
#module load fftw/3.3.4-gcc-8.4-openmpi-1.8-double
#module load fftw/3.3.4-gcc-8.4-openmpi-1.8-single
module load fftw/3.3.7-gcc-7.2-openmpi-3.0-single
module load fftw/3.3.7-gcc-7.2-openmpi-3.0-double

export CXX=mpicxx
export CC=mpicc
export FC=mpif90

export PYTHON_EXECUTABLE=python3

export LUSTREDIR=/cfs/klemming/nobackup/j/jnattila
export RUNKODIR=$LUSTREDIR/runko

export PYTHONPATH=$PYTHONPATH:$RUNKODIR
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/python
export PYTHONPATH=$PYTHONPATH:$RUNKODIR/analysis

export PYTHONDONTWRITEBYTECODE=true
export HDF5_USE_FILE_LOCKING=FALSE


# Tegner does not load all modules automatically (why?!) so we do this by hand
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pdc/vol/fftw/3.3.4/gcc/8.4/openmpi/1.8/single/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pdc/vol/fftw/3.3.4/gcc/8.4/openmpi/1.8/double/lib

#export PATH=$PATH:/pdc/vol/fftw/3.3.4/gcc/8.4/openmpi/1.8/single/include
#export PATH=$PATH:/pdc/vol/fftw/3.3.4/gcc/8.4/openmpi/1.8/double/include
