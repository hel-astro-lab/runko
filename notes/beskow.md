craype-hugepages512M


PrgEnv-intel
cdt/17.10
fftw
python/3.6.4
cray-hdf5/1.10.1.1
cmake/3.9.6
anaconda/py36/4.3
cray-mpich/7.6.3 cray-mpich/7.7.0
git

https://support.hdfgroup.org/ftp/HDF5/releases/ReleaseFiles/hdf5-1.10.1-RELEASE.txt


export HDF5_USE_FILE_LOCKING=FALSE


cdt/18.03
intel/18.0.0.128
module switch PrgEnv-cray/5.2.82 PrgEnv-intel
fftw/3.3.4.5
cray-hdf5/1.10.1.1
anaconda/py37/5.3
module switch gcc/4.8.1 gcc/7.3.0
cmake/3.12.1
craype-hugepages32M
git

source activate custom
conda activate mpi4py


pip install --user mpi4py




#setup miniconda environment
eval "$(/cfs/klemming/nobackup/j/jnattila/miniconda3/bin/conda shell.bash hook)"

conda install numpy
conda install mpi4py




anaconda/py36/4.3

anaconda/py36/4.3
source activate custom




#hack to install stuff on miniconda
pip install --prefix=/cfs/klemming/nobackup/j/jnattila/miniconda3 numpy
pip install --prefix=/cfs/klemming/nobackup/j/jnattila/miniconda3 mpi4py
