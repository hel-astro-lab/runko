#!/usr/bin/env bash

# Runko installation script for the HILE cluster (CPU build via preset).
# Run this on the login node (requires internet for downloads).

set -ve

#--------------------------------------------------
# Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

RUNKO_PATH=$(pwd)

#--------------------------------------------------
# Initialize git submodules and download rocThrust headers.
# rocThrust is needed for both GPU and CPU backends (tyvi uses thrust API).
# GPU builds find rocthrust via the rocm module; CPU builds use a minimal
# cmake config (cmake/rocthrust-cpu/) that points to these downloaded headers.
git submodule update --init --recursive

ROCM_LIBS_DIR="$RUNKO_PATH/external/tyvi/rocm-libraries"
if [ ! -d "$ROCM_LIBS_DIR/projects/rocthrust/thrust" ]; then
    cd "$RUNKO_PATH/external/tyvi"
    git clone --no-checkout --depth=1 --filter=tree:0 https://github.com/ROCm/rocm-libraries.git
    cd rocm-libraries
    git sparse-checkout init --cone
    git sparse-checkout set projects/rocthrust
    git checkout develop
    cd "$RUNKO_PATH"
fi

#--------------------------------------------------
# Pre-download mdspan (compute nodes have no internet; CPM needs this at configure time)
PRE_DOWNLOADED_MDSPAN_PATH="$HOME/mdspan"
if [ ! -d "$PRE_DOWNLOADED_MDSPAN_PATH" ]; then
    git clone https://github.com/kokkos/mdspan.git "$PRE_DOWNLOADED_MDSPAN_PATH"
fi

cd "$RUNKO_PATH"

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
# Create Python virtual environment for runko GPU builds
python -m venv venv/runko-gpu

# Activate the virtual environment
source venv/runko-gpu/bin/activate

# Install necessary Python dependencies
pip3 install h5py scipy matplotlib numpy

# Build mpi4py against Cray MPICH (PrgEnv-cray).
# --force-reinstall ensures the venv gets its own copy instead of using the
# system cray-python mpi4py, which is built against GNU MPICH.
MPI4PY_BUILD_MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

#--------------------------------------------------
# Patch the activate script with module loads, PYTHONPATH, and build helpers
P1="$RUNKO_PATH/"
P2="$RUNKO_PATH/external/corgi/lib"

cat >> venv/runko-gpu/bin/activate << EOL

# --- Runko environment (HILE GPU) ---
# Usage: source venv/runko-gpu/bin/activate

module purge
module load PrgEnv-cray
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load craype-x86-milan
module load cray-mpich
module load craype-network-ofi
module load rocm
module load libfabric
module load cray-python

export LD_LIBRARY_PATH=/opt/cray/pe/cce/18.0.1/cce/x86_64/lib:\$LD_LIBRARY_PATH

# GPU-aware MPI: preload GTL so it is available before mpi4py calls MPI_Init
export LD_PRELOAD=\${CRAY_MPICH_ROOTDIR}/gtl/lib/libmpi_gtl_hsa.so

# PYTHONPATH: repo root (runko_cpp_bindings.so) + corgi lib (pycorgi.so)
export PYTHONPATH="\$PYTHONPATH:${P1}:${P2}"
EOL


#--------------------------------------------------
# Next you need to build the code yourself!
#
# Get a compute node allocation:
#   srun -C g --gpus-per-node=1 --ntasks=1 --cpus-per-gpu=8 --mem=8G --time=03:00:00  --partition=hile --pty bash
#
# Activate the environment:
#   source venv/runko-gpu/bin/activate
#
# Configure and build using the preset:
#   cmake --preset hile-gpu-release -DCPM_mdspan_SOURCE=$HOME/mdspan
#   cmake --build hile-gpu-release -j8
