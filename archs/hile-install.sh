#!/usr/bin/env bash

# Runko installation script for on the HILE cluster.

set -v

#-------------------------------------------------- 
# Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

mkdir build

#-------------------------------------------------- 
# Pre-download mdspan (because there is no internet connection on compute node)
PRE_DOWNLOADED_MDSPAN_PATH="$HOME/mdspan"
git clone https://github.com/kokkos/mdspan.git $PRE_DOWNLOADED_MDSPAN_PATH



#-------------------------------------------------- 
# Update the PYTHONPATH environment variable with required runko
#  and corgi-related paths in order to make our venv runko-aware::
RUNKO_PATH=$(pwd)
P1="$RUNKO_PATH/"
P2="$RUNKO_PATH/external/corgi/lib"
export PYTHONPATH="$PYTHONPATH:$P1:$P2"

#-------------------------------------------------- 
# Load standard prerequisite modules for runko:
module load PrgEnv-cray
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load craype-x86-milan
module load cray-mpich 
module load craype-network-ofi
module load cray-python

#-------------------------------------------------- 
# Create a Python virtual environment specially for runko,
# located within the runko repository:
VENV_PATH="${RUNKO_PATH}/venvs"
mkdir $VENV_PATH
cd $VENV_PATH
python -m venv runko-gpu
cd $RUNKO_PATH

# Load this runko virtual environment:
source ${VENV_PATH}/runko-gpu/bin/activate

# Install necessary Python dependencies
pip3 install h5py scipy matplotlib numpy

# Build MPI4PY:
MPI4PY_BUILD_MPICC="cc -shared" python -m pip install --no-binary=mpi4py mpi4py


#-------------------------------------------------- 
# Next, we create a handy file which loads the necessary runko modules
# whenever called with "source archs/runko-load-env":
cat > runko-venv/bin/activate << EOL
# Tool to load runko modules
# Usage: "source runko-venv/bin/activate"
# Load standard prerequisite modules for runko:

module purge
module load PrgEnv-cray
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load craype-x86-milan
module load cray-mpich 
module load craype-network-ofi
module load rocm
module load libfabric
module load cray-python # not necessary as we are using a python virtual environment already
# Update the PYTHONPATH environment variable with runko path data
export PYTHONPATH="\$PYTHONPATH:${P1}:${P2}"
export CMAKE="/appl/lumi/SW/LUMI-24.03/common/EB/buildtools/24.03/bin/cmake -DCPM_mdspan_SOURCE=${PRE_DOWNLOADED_MDSPAN_PATH} -DCMAKE_CXX_COMPILER=CC"
EOL


#-------------------------------------------------- 
# Next you need to build the code yourself!
#
# Get gpu node allocation:
#   	srun -G1 -w hile-g02 --mem=8G --cpus-per-gpu=8 --time=03:00:00 --pty bash
#
# Load new modules:
# 	source archs/hile-load-runko-env
#
# Run cmake:
# 	cd build
# 	$CMAKE -DCMAKE_BUILD_TYPE=Release ..
#
# Finally, build:
# 	make -j8
