#!/usr/bin/env bash

# Runko installation script for a local macOS (Apple Silicon / M-series) machine.

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
# Create Python virtual environment for runko CPU builds
python3 -m venv venv/runko-v4

# Activate the virtual environment
source venv/runko-v4/bin/activate

# Install necessary Python dependencies
pip install --upgrade pip
pip install h5py scipy matplotlib numpy

# Build mpi4py against the Homebrew MPI that `mpicc` resolves to.
# --force-reinstall --no-cache-dir rebuilds mpi4py cleanly when the script
# is re-run (e.g. after a Homebrew MPI upgrade).
MPICC=mpicc pip install --no-cache-dir --no-binary=mpi4py mpi4py

#--------------------------------------------------
# Patch the activate script with PYTHONPATH so pycorgi / runko_cpp_bindings
P0="$RUNKO_PATH"
P1="$RUNKO_PATH/lib"
P2="$RUNKO_PATH/external/corgi/lib"

cat >> venv/runko-v4/bin/activate << EOL

# --- Runko environment (macOS CPU) ---
export PYTHONPATH="\$PYTHONPATH:${P0}:${P1}:${P2}"
EOL

#--------------------------------------------------
# Configure and build with the macOS vectorization-analysis preset.
export CC=mpicc
export CXX=mpic++

# Also export PYTHONPATH in the current shell so the check-runko target
# invoked by `make` can import pycorgi / pyrunko. The activate-script
# patch above only takes effect on *future* venv activations.
export PYTHONPATH="$PYTHONPATH:${P0}:${P1}:${P2}"

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPython_EXECUTABLE=$(which python) ..
make -j4


# source venv/runko-v4/bin/activate
