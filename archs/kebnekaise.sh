module load intel/2018b
module load HDF5
module load CMake
module load Eigen
module load Python/3.6.6

#intel compilers (icc is mpiicc and icpc is mpiicpc)
export CXX=mpiicpc 
export CC=mpiicc

#gnu compilers
#export CXX=mpic++
#export CC=mpicc

# python for CMAKE (susceptible for change)
export PYTHON_EXECUTABLE=/hpc2n/eb/software/MPI/intel/2018.3.222-GCC-7.3.0-2.30/impi/2018.3.222/Python/3.6.6/bin/python

# pybind installation
export pybind11_DIR=/home/n/natj/Public/libs/share/cmake/pybind11
export PYBIND11_DIR=/home/n/natj/Public/libs/share/cmake/pybind11

export EIGEN3_DIR=/home/n/natj/Public/libs/share/eigen3/cmake

# fmt installation fails wih intel compiler (14/17 compiler bug)
export FMT_INCLUDE_DIR=/home/n/natj/Public/fmt/include

export PYTHONDONTWRITEBYTECODE=true
export HDF5_USE_FILE_LOCKING=FALSE

export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/runko
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/runko/lib
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/runko/corgi/lib
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/runko/python
export PYTHONPATH=$PYTHONPATH:/pfs/nobackup/home/n/natj/runko/analysis

