#!/usr/bin/env bash

# Runko installation script for the LUMI supercomputer (CPU build).
# Mirrors the structure of macos-install.sh; uses the Cray PrgEnv + cray-python
# stack instead of Homebrew MPI.

set -ve

#--------------------------------------------------
# Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

RUNKO_PATH=$(pwd)

#--------------------------------------------------
# Initialize git submodules and download rocThrust headers.
# rocThrust is needed for the CPU backend (tyvi uses the thrust API);
# CPU builds resolve it via the minimal cmake config (cmake/rocthrust-cpu/)
# that points to these downloaded headers.
git submodule update --init --recursive
cd "$RUNKO_PATH"

#--------------------------------------------------
# Load LUMI prerequisite modules for a CPU build.
# NOTE: rocm + craype-accel-amd-gfx90a are intentionally NOT loaded — with them
# the Cray CC wrapper would compile everything as HIP/GPU code.
module load LUMI
module load partition/C
module load PrgEnv-cray
module load cray-hdf5
module load cray-mpich craype-network-ofi
module load buildtools
module load cray-python

#--------------------------------------------------
# Create Python virtual environment for runko CPU builds
python -m venv venv/runko-v4

# Activate the virtual environment
source venv/runko-v4/bin/activate

# Install necessary Python dependencies
pip install --upgrade pip
pip install h5py scipy matplotlib numpy

# Build mpi4py against Cray MPICH (PrgEnv-cray).
# --force-reinstall ensures the venv gets its own copy instead of shadowing
# the cray-python mpi4py (which is built against GNU MPICH).
MPI4PY_BUILD_MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

#--------------------------------------------------
# Patch the activate script so future `source venv/runko-v4/bin/activate`
# calls re-load the module set and set PYTHONPATH for pycorgi / pyrunko /
# pytools.
#   P0 = repo root (for the `pytools` python package)
#   P1 = pyrunko .so location
#   P2 = pycorgi .so location
P0="$RUNKO_PATH"
P1="$RUNKO_PATH/lib"
P2="$RUNKO_PATH/external/corgi/lib"

cat >> venv/runko-v4/bin/activate << EOL

# --- Runko environment (LUMI CPU) ---
module load LUMI
module load partition/C
module load PrgEnv-cray
module load cray-hdf5
module load cray-mpich craype-network-ofi
module load buildtools

export PYTHONPATH="\$PYTHONPATH:${P0}:${P1}:${P2}"
EOL

#--------------------------------------------------
# Configure and build with the Cray compiler wrappers.
export CC=cc
export CXX=CC

# Also export PYTHONPATH in the current shell so the check-runko target
# invoked by `make` can import pycorgi / pyrunko / pytools. The activate-script
# patch above only takes effect on *future* venv activations.
export PYTHONPATH="$PYTHONPATH:${P0}:${P1}:${P2}"

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPython_EXECUTABLE=$(which python) ..
make -j18


# source venv/runko-v4/bin/activate
